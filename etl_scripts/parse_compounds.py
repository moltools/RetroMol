"""ETL for adding retrieval objects and primary sequences to compounds using RetroMol stream.

Default behavior:
    Process only compounds that do not yet have a RetrievalObject.

Force behavior:
    Recalculate all selected compounds by deleting existing RetrievalObject rows.
    PrimarySequence rows are deleted automatically by cascade.

RetroMol behavior:
    Compounds are streamed through run_retromol_stream(), which can parallelize
    parsing over multiple workers.
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from typing import Any, Iterator

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import insert as pg_insert
from sqlalchemy.exc import SQLAlchemyError
from tqdm import tqdm
from rdkit import RDLogger

from retromol.io.streaming import ResultEvent, run_retromol_stream
from retromol.model.result import Result
from retromol.model.rules import RuleSet

from retromol_database.db.engine import SessionLocal
from retromol_database.db.models import (
    BIOSYNTHETIC_FP_SIZE,
    Compound,
    PrimarySequence,
    PrimarySequenceMonomer,
    RetrievalObject,
)


RDLogger.DisableLog("rdApp.*")


log = logging.getLogger(__name__)


@dataclass(frozen=True)
class CompoundInput:
    """Minimal compound input needed for RetroMol parsing."""

    id: int
    inchikey: str
    smiles: str


@dataclass(frozen=True)
class CalculatedPrimarySequence:
    """Calculated primary sequence for one compound."""

    tokens: list[str]
    source: str = "retromol"
    label: str | None = None

    # Sequence-specific fingerprint.
    biosynthetic_fp: list[float] | None = None
    biosynthetic_fp_model: str | None = None

    score: float | None = None
    extra: dict[str, Any] | None = None


@dataclass(frozen=True)
class CalculatedCompoundRetrieval:
    """Calculated retrieval data for one compound."""

    # Unordered bag of all identified monomers.
    bag_of_monomers: list[str]
    bag_of_monomers_model: str | None

    coverage: float | None
    primary_sequences: list[CalculatedPrimarySequence]

    extra: dict[str, Any] | None = None


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--batch-size", type=int, default=500)
    parser.add_argument("--limit", type=int, default=None)

    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of RetroMol worker processes.",
    )
    parser.add_argument(
        "--retromol-batch-size",
        type=int,
        default=2_000,
        help="Number of compounds buffered by run_retromol_stream.",
    )
    parser.add_argument(
        "--pool-chunksize",
        type=int,
        default=50,
        help="Multiprocessing imap_unordered chunksize for RetroMol workers.",
    )
    parser.add_argument(
        "--maxtasksperchild",
        type=int,
        default=2_000,
        help="Restart RetroMol workers after this many tasks.",
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help="Delete existing retrieval objects and primary sequences before recalculating.",
    )

    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Run RetroMol and parsing but do not write to the database.",
    )

    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


# =============================================================================
# Result parsing
# =============================================================================


def parse_retromol_result(
    result: Result,
    compound: CompoundInput,
) -> CalculatedCompoundRetrieval:
    """Convert one RetroMol Result into database-ready retrieval output.

    Example skeleton:

        coverage = result.calculate_coverage()

        bag_of_monomers = [...]
        primary_sequences = [
            CalculatedPrimarySequence(
                tokens=[...],
                source="retromol",
                label="best_parse",
                biosynthetic_fp=[...],
                biosynthetic_fp_model="retromol_sequence_fp_v1",
                score=None,
                extra={...},
            )
        ]

        return CalculatedCompoundRetrieval(
            bag_of_monomers=bag_of_monomers,
            bag_of_monomers_model="retromol_bag_v1",
            coverage=coverage,
            primary_sequences=primary_sequences,
            extra={"retromol": result.to_dict()},
        )
    """
    print(result)
    print(compound)

    raise NotImplementedError


def validate_calculated_retrieval(result: CalculatedCompoundRetrieval) -> CalculatedCompoundRetrieval:
    """Validate and normalize calculated retrieval output."""

    for seq in result.primary_sequences:
        if seq.biosynthetic_fp is not None and len(seq.biosynthetic_fp) != BIOSYNTHETIC_FP_SIZE:
            raise ValueError(
                f"primary sequence biosynthetic_fp has length {len(seq.biosynthetic_fp)}, "
                f"expected {BIOSYNTHETIC_FP_SIZE}"
            )

    return CalculatedCompoundRetrieval(
        bag_of_monomers=[str(x) for x in result.bag_of_monomers],
        bag_of_monomers_model=result.bag_of_monomers_model,
        coverage=result.coverage,
        primary_sequences=result.primary_sequences,
        extra=result.extra,
    )


# =============================================================================
# RetroMol stream helpers
# =============================================================================


def compounds_to_retromol_rows(compounds: list[CompoundInput]) -> Iterator[dict[str, Any]]:
    """Convert compounds into rows accepted by run_retromol_stream.

    run_retromol_stream passes the full row dict as Submission.props.
    That means compound_id comes back inside result.submission.props.
    """

    for compound in compounds:
        yield {
            "compound_id": compound.id,
            "inchikey": compound.inchikey,
            "smiles": compound.smiles,
        }


def run_retromol_for_compounds(
    compounds: list[CompoundInput],
    ruleset: RuleSet,
    workers: int,
    retromol_batch_size: int,
    pool_chunksize: int,
    maxtasksperchild: int,
) -> Iterator[tuple[CompoundInput, ResultEvent]]:
    """Run RetroMol stream and yield events mapped back to CompoundInput."""

    compound_by_id = {compound.id: compound for compound in compounds}

    for evt in run_retromol_stream(
        ruleset=ruleset,
        row_iter=compounds_to_retromol_rows(compounds),
        smiles_col="smiles",
        workers=workers,
        batch_size=retromol_batch_size,
        pool_chunksize=pool_chunksize,
        maxtasksperchild=maxtasksperchild,
    ):
        compound_id: int | None = None

        if evt.result is not None:
            try:
                result = Result.from_dict(evt.result)
                props = result.submission.props or {}
                compound_id = int(props["compound_id"])
            except Exception as e:
                log.warning("Could not map RetroMol result back to compound_id: %s", e)
                continue

        # If RetroMol errors before returning a serialized Result, the current
        # ResultEvent does not include props. In that case we cannot reliably map
        # the error to a compound unless run_retromol_stream is changed to return
        # props with errors too.
        if compound_id is None:
            log.warning("RetroMol failed with unmapped error: %s", evt.error)
            continue

        compound = compound_by_id.get(compound_id)
        if compound is None:
            log.warning("RetroMol returned unknown compound_id=%s", compound_id)
            continue

        yield compound, evt


# =============================================================================
# Database helpers
# =============================================================================


def iter_compound_id_batches(
    batch_size: int,
    force: bool,
    limit: int | None,
) -> Iterator[list[int]]:
    """Yield compound IDs to process.

    If force=False:
        only compounds without a RetrievalObject are selected.

    If force=True:
        all compounds are selected.
    """

    yielded = 0
    last_id = 0

    while True:
        with SessionLocal() as session:
            stmt = sa.select(Compound.id).where(Compound.id > last_id)

            if not force:
                stmt = (
                    stmt.outerjoin(
                        RetrievalObject,
                        RetrievalObject.compound_id == Compound.id,
                    )
                    .where(RetrievalObject.id.is_(None))
                )

            stmt = stmt.order_by(Compound.id).limit(batch_size)

            if limit is not None:
                remaining = limit - yielded
                if remaining <= 0:
                    break
                stmt = stmt.limit(min(batch_size, remaining))

            ids = list(session.execute(stmt).scalars().all())

        if not ids:
            break

        yielded += len(ids)
        last_id = ids[-1]

        yield ids


def load_compound_inputs_by_ids(session, compound_ids: list[int]) -> list[CompoundInput]:
    """Load minimal compound data for a batch of IDs."""

    rows = session.execute(
        sa.select(
            Compound.id,
            Compound.inchikey,
            Compound.smiles,
        )
        .where(Compound.id.in_(compound_ids))
        .order_by(Compound.id)
    ).all()

    return [
        CompoundInput(
            id=compound_id,
            inchikey=inchikey,
            smiles=smiles,
        )
        for compound_id, inchikey, smiles in rows
    ]


def delete_existing_outputs(session, compound_ids: list[int]) -> None:
    """Delete existing retrieval objects for compounds.

    PrimarySequence and PrimarySequenceMonomer rows are deleted by cascade.
    """

    if not compound_ids:
        return

    session.execute(
        sa.delete(RetrievalObject).where(
            RetrievalObject.compound_id.in_(compound_ids)
        )
    )
    session.flush()


def upsert_retrieval_object(
    session,
    compound_id: int,
    result: CalculatedCompoundRetrieval,
) -> int:
    """Insert/update retrieval object and return retrieval_object.id."""

    row = {
        "compound_id": compound_id,
        "gene_cluster_id": None,
        "bag_of_monomers": result.bag_of_monomers,
        "bag_of_monomers_model": result.bag_of_monomers_model,
        "coverage": result.coverage,
        "extra": result.extra,
    }

    stmt = pg_insert(RetrievalObject).values(row)

    session.execute(
        stmt.on_conflict_do_update(
            index_elements=[RetrievalObject.compound_id],
            set_={
                "bag_of_monomers": stmt.excluded.bag_of_monomers,
                "bag_of_monomers_model": stmt.excluded.bag_of_monomers_model,
                "coverage": stmt.excluded.coverage,
                "extra": stmt.excluded.extra,
            },
        )
    )
    session.flush()

    retrieval_object_id = session.execute(
        sa.select(RetrievalObject.id).where(
            RetrievalObject.compound_id == compound_id
        )
    ).scalar_one()

    return retrieval_object_id


def insert_primary_sequences(
    session,
    retrieval_object_id: int,
    sequences: list[CalculatedPrimarySequence],
) -> None:
    """Insert primary sequences and their monomer-position rows."""

    for seq in sequences:
        tokens = [str(token) for token in seq.tokens]

        primary_sequence = PrimarySequence(
            retrieval_object_id=retrieval_object_id,
            source=seq.source,
            label=seq.label,
            length=len(tokens),
            tokens=tokens,
            biosynthetic_fp=seq.biosynthetic_fp,
            biosynthetic_fp_model=seq.biosynthetic_fp_model,
            score=seq.score,
            extra=seq.extra,
        )

        session.add(primary_sequence)
        session.flush()

        monomer_rows = [
            {
                "sequence_id": primary_sequence.id,
                "position": idx,
                "token": token,
                "token_class": None,
                "raw_token": token,
                "confidence": None,
                "extra": None,
            }
            for idx, token in enumerate(tokens)
        ]

        if monomer_rows:
            session.execute(pg_insert(PrimarySequenceMonomer).values(monomer_rows))


def write_compound_outputs(
    session,
    compound_id: int,
    result: CalculatedCompoundRetrieval,
) -> None:
    """Write calculated retrieval object and primary sequences for one compound."""

    retrieval_object_id = upsert_retrieval_object(session, compound_id, result)

    # If retrieval object already existed, replace its previous sequences.
    session.execute(
        sa.delete(PrimarySequence).where(
            PrimarySequence.retrieval_object_id == retrieval_object_id
        )
    )
    session.flush()

    insert_primary_sequences(
        session=session,
        retrieval_object_id=retrieval_object_id,
        sequences=result.primary_sequences,
    )


# =============================================================================
# Main ETL
# =============================================================================


def run(
    batch_size: int,
    workers: int,
    retromol_batch_size: int,
    pool_chunksize: int,
    maxtasksperchild: int,
    force: bool,
    limit: int | None,
    dry_run: bool,
) -> None:
    """Run ETL."""

    ruleset = RuleSet.load_default()

    total_seen = 0
    total_success = 0
    total_failed = 0
    total_retromol_errors = 0
    total_sequences = 0
    total_monomers = 0

    batches = iter_compound_id_batches(
        batch_size=batch_size,
        force=force,
        limit=limit,
    )

    for compound_ids in tqdm(batches, desc="Processing compound batches"):
        total_seen += len(compound_ids)

        try:
            with SessionLocal() as session:
                compounds = load_compound_inputs_by_ids(session, compound_ids)
        except SQLAlchemyError as e:
            total_failed += len(compound_ids)
            log.error(
                "Database error loading batch starting compound_id=%s: %s",
                compound_ids[0],
                e,
            )
            continue

        calculated: list[tuple[int, CalculatedCompoundRetrieval]] = []

        for compound, evt in run_retromol_for_compounds(
            compounds=compounds,
            ruleset=ruleset,
            workers=workers,
            retromol_batch_size=retromol_batch_size,
            pool_chunksize=pool_chunksize,
            maxtasksperchild=maxtasksperchild,
        ):
            if evt.error is not None:
                total_retromol_errors += 1
                total_failed += 1
                log.warning(
                    "RetroMol failed for compound id=%s inchikey=%s: %s",
                    compound.id,
                    compound.inchikey,
                    evt.error,
                )
                continue

            if evt.result is None:
                total_failed += 1
                log.warning(
                    "RetroMol returned empty result for compound id=%s inchikey=%s",
                    compound.id,
                    compound.inchikey,
                )
                continue

            try:
                retromol_result = Result.from_dict(evt.result)

                result = validate_calculated_retrieval(
                    parse_retromol_result(
                        result=retromol_result,
                        compound=compound,
                    )
                )

                calculated.append((compound.id, result))

                total_success += 1
                total_sequences += len(result.primary_sequences)
                total_monomers += len(result.bag_of_monomers)

            except NotImplementedError:
                raise

            except Exception as e:
                total_failed += 1
                log.warning(
                    "Failed to parse RetroMol output for compound id=%s inchikey=%s: %s",
                    compound.id,
                    compound.inchikey,
                    e,
                )

        if dry_run:
            continue

        try:
            with SessionLocal() as session:
                if force:
                    delete_existing_outputs(session, compound_ids)

                for compound_id, result in calculated:
                    write_compound_outputs(
                        session=session,
                        compound_id=compound_id,
                        result=result,
                    )

                session.commit()

        except SQLAlchemyError as e:
            total_failed += len(calculated)
            log.error(
                "Database error writing batch starting compound_id=%s: %s",
                compound_ids[0],
                e,
            )

    log.info(
        "Finished compound retrieval ETL. "
        "seen=%s success=%s failed=%s retromol_errors=%s "
        "primary_sequences=%s bag_monomers=%s force=%s dry_run=%s workers=%s",
        total_seen,
        total_success,
        total_failed,
        total_retromol_errors,
        total_sequences,
        total_monomers,
        force,
        dry_run,
        workers,
    )


def main() -> None:
    args = cli()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    run(
        batch_size=args.batch_size,
        workers=args.workers,
        retromol_batch_size=args.retromol_batch_size,
        pool_chunksize=args.pool_chunksize,
        maxtasksperchild=args.maxtasksperchild,
        force=args.force,
        limit=args.limit,
        dry_run=args.dry_run,
    )


if __name__ == "__main__":
    main()