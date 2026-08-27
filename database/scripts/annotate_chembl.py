"""Step 10: annotate every compound entry's bioactivity via ChEMBL.

Runs after compound loading (like annotate_npclassifier.py), so every distinct
molecule -- regardless of which source(s) it came from -- is looked up exactly once.
BGC entries are skipped: bioactivity here is a property of the compound structure
(looked up by its InChIKey entry id), not of the producing organism/cluster.

Unlike NPClassifier, this queries a local ChEMBL SQLite release (see chembl.py) rather
than a rate-limited API -- no pacing or caching needed, lookups are just local reads.
"""

import argparse
import logging
from pathlib import Path

from tqdm import tqdm

from chembl import ChemblDB
from retromol_database.duckdb import RetroMolDuckDB

log = logging.getLogger(__name__)


def run(db_path: str | Path, chembl_sqlite_path: str | Path, log_every: int = 500) -> None:
    processed = 0
    matched = 0
    annotated = 0

    db = RetroMolDuckDB.open(db_path)
    try:
        total = db.count_entries_by_type(["compound"])
        with ChemblDB(chembl_sqlite_path) as chembl:
            with tqdm(total=total, desc="annotate_chembl", unit="cmpd") as pbar:
                for entry in db.iter_entries():
                    if entry.type != "compound":
                        continue

                    processed += 1
                    result = chembl.bioactivity_for_inchikey(entry.id)
                    if result is not None:
                        matched += 1
                        if result.max_phase_label:
                            db.add_bioactivity_annotation(
                                entry.id,
                                level="chembl_max_phase",
                                label=result.max_phase_label,
                                external_id=result.chembl_id,
                            )
                            annotated += 1
                        for atc in result.atc_categories:
                            db.add_bioactivity_annotation(
                                entry.id,
                                level="chembl_atc",
                                label=atc.level1_description,
                                external_id=atc.level5,
                            )
                            annotated += 1

                    pbar.update(1)
                    pbar.set_postfix(matched=matched, annotated=annotated)

                    if log_every > 0 and processed % log_every == 0:
                        log.info(
                            "annotate_chembl: processed=%d matched=%d annotated=%d", processed, matched, annotated
                        )
    finally:
        db.close()

    log.info("annotate_chembl: processed=%d matched=%d annotated=%d", processed, matched, annotated)


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--chembl-sqlite-path", required=True)
    ap.add_argument("--log-every", type=int, default=500)
    args = ap.parse_args()

    run(db_path=args.db_path, chembl_sqlite_path=args.chembl_sqlite_path, log_every=args.log_every)


if __name__ == "__main__":
    main()
