"""ETL for loading NPAtlas compounds, compound names, and compound database references."""

from __future__ import annotations

import argparse
import json
import logging
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterator, TypeVar

import sqlalchemy as sa
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from sqlalchemy.dialects.postgresql import insert as pg_insert
from sqlalchemy.exc import SQLAlchemyError
from tqdm import tqdm

from retromol_database.db.engine import SessionLocal
from retromol_database.db.models import (
    Compound,
    CompoundName,
    DatabaseReference,
    compound_reference,
)
from retromol_database.utils.download import download_and_prepare


RDLogger.DisableLog("rdApp.*")

NPATLAS_URL = r"https://www.npatlas.org/static/downloads/NPAtlas_download.json"

MORGAN_RADIUS = 2
MORGAN_N_BITS = 2048

log = logging.getLogger(__name__)

T = TypeVar("T")


@dataclass(frozen=True)
class CompoundProps:
    """Basic compound properties calculated from RDKit."""

    formula: str | None
    mol_weight: float | None
    morgan_fp: str | None
    c_atom_count: int | None
    h_atom_count: int | None
    n_atom_count: int | None
    o_atom_count: int | None
    p_atom_count: int | None
    s_atom_count: int | None
    f_atom_count: int | None
    cl_atom_count: int | None
    br_atom_count: int | None
    i_atom_count: int | None


@dataclass(frozen=True)
class NPAtlasCompoundJob:
    """Single NPAtlas compound parsed from the NPAtlas JSON."""

    inchikey: str
    smiles: str
    npaid: str
    original_name: str | None
    synonyms: list[str]


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=1_000)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--log-level", default="INFO")
    return parser.parse_args()


def iter_chunks(items: list[T], chunk_size: int) -> Iterator[list[T]]:
    """Yield chunks from a list."""

    for i in range(0, len(items), chunk_size):
        yield items[i : i + chunk_size]


def normalize_text(value: Any, max_len: int | None = None) -> str | None:
    """Normalize optional text from NPAtlas."""

    if value is None:
        return None

    if isinstance(value, dict):
        value = value.get("name")

    value = str(value).strip()
    if not value:
        return None

    if max_len is not None:
        value = value[:max_len]

    return value


def normalize_synonyms(value: Any) -> list[str]:
    """Normalize NPAtlas synonyms into a clean unique list."""

    if not value:
        return []

    if not isinstance(value, list):
        value = [value]

    out: list[str] = []
    seen: set[str] = set()

    for item in value:
        name = normalize_text(item)
        if name is None:
            continue

        key = name.lower()
        if key in seen:
            continue

        seen.add(key)
        out.append(name)

    return out


def parse_formula_counts(formula: str | None) -> dict[str, int | None]:
    """Parse atom counts from a formula string.

    Example:
        C10H15N -> {"C": 10, "H": 15, "N": 1}

    This is enough for the tracked elements in the schema.
    """

    elements = {
        "C": 0,
        "H": 0,
        "N": 0,
        "O": 0,
        "P": 0,
        "S": 0,
        "F": 0,
        "Cl": 0,
        "Br": 0,
        "I": 0,
    }

    if not formula:
        return {k: None for k in elements}

    for element, count in re.findall(r"([A-Z][a-z]?)(\d*)", formula):
        if element not in elements:
            continue

        elements[element] += int(count) if count else 1

    return elements


def calculate_compound_props(smiles: str) -> CompoundProps:
    """Calculate basic compound properties and Morgan fingerprint."""

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")

    formula = rdMolDescriptors.CalcMolFormula(mol)
    counts = parse_formula_counts(formula)

    fp = AllChem.GetMorganFingerprintAsBitVect(
        mol,
        radius=MORGAN_RADIUS,
        nBits=MORGAN_N_BITS,
    )
    morgan_fp = DataStructs.BitVectToText(fp)

    return CompoundProps(
        formula=formula,
        mol_weight=float(Descriptors.MolWt(mol)),
        morgan_fp=morgan_fp,
        c_atom_count=counts["C"],
        h_atom_count=counts["H"],
        n_atom_count=counts["N"],
        o_atom_count=counts["O"],
        p_atom_count=counts["P"],
        s_atom_count=counts["S"],
        f_atom_count=counts["F"],
        cl_atom_count=counts["Cl"],
        br_atom_count=counts["Br"],
        i_atom_count=counts["I"],
    )


def parse_npatlas_jobs(data: list[dict[str, Any]]) -> list[NPAtlasCompoundJob]:
    """Parse NPAtlas JSON records into jobs."""

    jobs: list[NPAtlasCompoundJob] = []
    seen_inchikeys: set[str] = set()

    for item in tqdm(data, desc="Parsing NPAtlas JSON"):
        inchikey = normalize_text(item.get("inchikey"))
        smiles = normalize_text(item.get("smiles"))
        npaid = normalize_text(item.get("npaid"))
        original_name = normalize_text(item.get("original_name"))
        synonyms = normalize_synonyms(item.get("synonyms", []))

        if inchikey is None or smiles is None or npaid is None:
            log.warning("Skipping malformed NPAtlas entry missing inchikey/smiles/npaid.")
            continue

        if inchikey in seen_inchikeys:
            log.warning("Skipping duplicate InChIKey within NPAtlas JSON: %s", inchikey)
            continue

        seen_inchikeys.add(inchikey)

        jobs.append(
            NPAtlasCompoundJob(
                inchikey=inchikey,
                smiles=smiles,
                npaid=npaid,
                original_name=original_name,
                synonyms=synonyms,
            )
        )

    return jobs


def compound_row_for_job(job: NPAtlasCompoundJob) -> dict[str, Any]:
    """Build a Compound insert row from an NPAtlas job."""

    props = calculate_compound_props(job.smiles)

    return {
        "inchikey": job.inchikey,
        "smiles": job.smiles,
        "canonical_name": job.original_name,
        "formula": props.formula,
        "mol_weight": props.mol_weight,
        "c_atom_count": props.c_atom_count,
        "h_atom_count": props.h_atom_count,
        "n_atom_count": props.n_atom_count,
        "o_atom_count": props.o_atom_count,
        "p_atom_count": props.p_atom_count,
        "s_atom_count": props.s_atom_count,
        "f_atom_count": props.f_atom_count,
        "cl_atom_count": props.cl_atom_count,
        "br_atom_count": props.br_atom_count,
        "i_atom_count": props.i_atom_count,
        "morgan_fp": props.morgan_fp,
        "morgan_fp_radius": MORGAN_RADIUS,
        "morgan_fp_n_bits": MORGAN_N_BITS,
        "parser_payload": None,
    }


def compound_name_rows_for_job(job: NPAtlasCompoundJob, compound_id: int) -> list[dict[str, Any]]:
    """Build CompoundName rows for a job."""

    rows: list[dict[str, Any]] = []
    seen: set[str] = set()

    def add_name(name: str | None, name_type: str) -> None:
        if name is None:
            return

        key = name.lower()
        if key in seen:
            return

        seen.add(key)
        rows.append(
            {
                "compound_id": compound_id,
                "name": name,
                "name_type": name_type,
                "source": "NPAtlas",
            }
        )

    add_name(job.original_name, "original")

    for synonym in job.synonyms:
        add_name(synonym, "synonym")

    return rows


def database_reference_row_for_job(job: NPAtlasCompoundJob) -> dict[str, Any]:
    """Build DatabaseReference row for the NPAtlas compound accession."""

    return {
        "database_name": "NPAtlas",
        "database_identifier": job.npaid,
        "name": job.original_name,
        "url": None,
        "extra": None,
    }


def upsert_compounds_and_get_ids(session, compound_rows: list[dict[str, Any]]) -> dict[str, int]:
    """Upsert compounds and return inchikey -> compound_id."""

    if not compound_rows:
        return {}

    stmt = pg_insert(Compound).values(compound_rows)

    update_cols = {
        "smiles": stmt.excluded.smiles,
        "canonical_name": stmt.excluded.canonical_name,
        "formula": stmt.excluded.formula,
        "mol_weight": stmt.excluded.mol_weight,
        "c_atom_count": stmt.excluded.c_atom_count,
        "h_atom_count": stmt.excluded.h_atom_count,
        "n_atom_count": stmt.excluded.n_atom_count,
        "o_atom_count": stmt.excluded.o_atom_count,
        "p_atom_count": stmt.excluded.p_atom_count,
        "s_atom_count": stmt.excluded.s_atom_count,
        "f_atom_count": stmt.excluded.f_atom_count,
        "cl_atom_count": stmt.excluded.cl_atom_count,
        "br_atom_count": stmt.excluded.br_atom_count,
        "i_atom_count": stmt.excluded.i_atom_count,
        "morgan_fp": stmt.excluded.morgan_fp,
        "morgan_fp_radius": stmt.excluded.morgan_fp_radius,
        "morgan_fp_n_bits": stmt.excluded.morgan_fp_n_bits,
        "parser_payload": stmt.excluded.parser_payload,
    }

    session.execute(
        stmt.on_conflict_do_update(
            index_elements=[Compound.inchikey],
            set_=update_cols,
        )
    )
    session.flush()

    inchikeys = [row["inchikey"] for row in compound_rows]
    rows = session.execute(
        sa.select(Compound.inchikey, Compound.id).where(Compound.inchikey.in_(inchikeys))
    ).all()

    return {inchikey: compound_id for inchikey, compound_id in rows}


def upsert_database_references_and_get_ids(
    session,
    reference_rows: list[dict[str, Any]],
) -> dict[tuple[str, str], int]:
    """Upsert database references and return (database_name, database_identifier) -> id."""

    if not reference_rows:
        return {}

    stmt = pg_insert(DatabaseReference).values(reference_rows)

    session.execute(
        stmt.on_conflict_do_update(
            constraint="ux_database_reference_database_identifier",
            set_={
                "name": stmt.excluded.name,
                "url": stmt.excluded.url,
                "extra": stmt.excluded.extra,
            },
        )
    )
    session.flush()

    keys = {
        (row["database_name"], row["database_identifier"])
        for row in reference_rows
    }

    rows = session.execute(
        sa.select(
            DatabaseReference.database_name,
            DatabaseReference.database_identifier,
            DatabaseReference.id,
        ).where(
            sa.tuple_(
                DatabaseReference.database_name,
                DatabaseReference.database_identifier,
            ).in_(list(keys))
        )
    ).all()

    return {
        (database_name, database_identifier): reference_id
        for database_name, database_identifier, reference_id in rows
    }


def insert_compound_names(session, rows: list[dict[str, Any]]) -> None:
    """Insert compound names, ignoring duplicates."""

    if not rows:
        return

    stmt = (
        pg_insert(CompoundName)
        .values(rows)
        .on_conflict_do_nothing(
            constraint="ux_compound_name_compound_name_source",
        )
    )
    session.execute(stmt)


def insert_compound_reference_links(session, rows: list[dict[str, Any]]) -> None:
    """Insert compound-reference links, ignoring duplicates."""

    if not rows:
        return

    stmt = (
        pg_insert(compound_reference)
        .values(rows)
        .on_conflict_do_nothing(
            index_elements=[
                compound_reference.c.compound_id,
                compound_reference.c.reference_id,
            ]
        )
    )
    session.execute(stmt)


def load_compounds_npatlas(
    tmp: Path,
    batch_size: int = 1_000,
    limit: int | None = None,
) -> None:
    """Load NPAtlas compounds, names, and database references."""

    tmp = tmp.expanduser()
    tmp.mkdir(parents=True, exist_ok=True)

    npatlas_path = download_and_prepare(url=NPATLAS_URL, base_cache=tmp)

    with open(npatlas_path, "r") as f_i:
        data = json.load(f_i)

    jobs = parse_npatlas_jobs(data)

    if limit is not None:
        jobs = jobs[:limit]

    log.info("Prepared %s NPAtlas compound jobs.", len(jobs))

    total_compounds = 0
    total_names = 0
    total_references = 0
    total_links = 0
    failed = 0

    for batch_jobs in tqdm(
            iter_chunks(jobs, batch_size),
            total=(len(jobs) + batch_size - 1) // batch_size,
            desc="Loading NPAtlas compound batches",
    ):
        compound_rows: list[dict[str, Any]] = []
        valid_jobs: list[NPAtlasCompoundJob] = []

        for job in batch_jobs:
            try:
                compound_rows.append(compound_row_for_job(job))
                valid_jobs.append(job)
            except Exception as e:
                failed += 1
                log.warning("Skipping NPAtlas compound npaid=%s inchikey=%s: %s", job.npaid, job.inchikey, e)

        if not compound_rows:
            continue

        try:
            with SessionLocal() as session:
                inchikey_to_compound_id = upsert_compounds_and_get_ids(session, compound_rows)

                reference_rows = [
                    database_reference_row_for_job(job)
                    for job in valid_jobs
                    if job.inchikey in inchikey_to_compound_id
                ]

                reference_key_to_id = upsert_database_references_and_get_ids(
                    session,
                    reference_rows,
                )

                name_rows: list[dict[str, Any]] = []
                compound_reference_rows: list[dict[str, Any]] = []

                for job in valid_jobs:
                    compound_id = inchikey_to_compound_id.get(job.inchikey)
                    if compound_id is None:
                        continue

                    name_rows.extend(compound_name_rows_for_job(job, compound_id))

                    reference_id = reference_key_to_id.get(("NPAtlas", job.npaid))
                    if reference_id is not None:
                        compound_reference_rows.append(
                            {
                                "compound_id": compound_id,
                                "reference_id": reference_id,
                            }
                        )

                insert_compound_names(session, name_rows)
                insert_compound_reference_links(session, compound_reference_rows)

                session.commit()

                total_compounds += len(inchikey_to_compound_id)
                total_names += len(name_rows)
                total_references += len(reference_rows)
                total_links += len(compound_reference_rows)

        except SQLAlchemyError as e:
            failed += len(compound_rows)
            log.error("Database error during NPAtlas batch flush: %s", e)

    log.info(
        "Finished loading NPAtlas compounds. "
        "compounds=%s, names=%s, references=%s, links=%s, failed=%s",
        total_compounds,
        total_names,
        total_references,
        total_links,
        failed,
    )


def main() -> None:
    args = cli()

    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    load_compounds_npatlas(
        tmp=args.tmp,
        batch_size=args.batch_size,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()