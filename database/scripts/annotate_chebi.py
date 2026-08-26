"""Step 11: annotate every compound entry's bioactivity via ChEBI's role ontology.

Runs after compound loading, so every distinct molecule -- regardless of which
source(s) it came from -- is looked up exactly once. BGC entries are skipped, same
reasoning as annotate_chembl.py.

Like annotate_chembl.py, this queries a local bulk release (see chebi.py) rather than
a rate-limited API -- no pacing or caching needed.
"""

import argparse
import logging
from pathlib import Path

from chebi import ChebiDB
from retromol_database.duckdb import RetroMolDuckDB

log = logging.getLogger(__name__)


def run(db_path: str | Path, chebi_dir: str | Path, log_every: int = 500) -> None:
    processed = 0
    matched = 0
    annotated = 0

    db = RetroMolDuckDB.open(db_path)
    try:
        chebi = ChebiDB.load(chebi_dir)
        for entry in db.iter_entries():
            if entry.type != "compound":
                continue

            processed += 1
            result = chebi.roles_for_inchikey(entry.id)
            if result is not None:
                matched += 1
                for role in result.biological_roles:
                    db.add_bioactivity_annotation(
                        entry.id, level="chebi_biological_role", label=role.label, external_id=role.chebi_accession
                    )
                    annotated += 1
                for role in result.chemical_roles:
                    db.add_bioactivity_annotation(
                        entry.id, level="chebi_chemical_role", label=role.label, external_id=role.chebi_accession
                    )
                    annotated += 1

            if log_every > 0 and processed % log_every == 0:
                log.info("annotate_chebi: processed=%d matched=%d annotated=%d", processed, matched, annotated)
    finally:
        db.close()

    log.info("annotate_chebi: processed=%d matched=%d annotated=%d", processed, matched, annotated)


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--chebi-dir", required=True, help="dir with ChEBI's compounds.tsv.gz/structures.tsv.gz/relation.tsv.gz")
    ap.add_argument("--log-every", type=int, default=500)
    args = ap.parse_args()

    run(db_path=args.db_path, chebi_dir=args.chebi_dir, log_every=args.log_every)


if __name__ == "__main__":
    main()
