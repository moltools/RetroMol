"""Step 1: create an empty RetroMol database."""

import argparse
from pathlib import Path

from retromol_database.duckdb import RetroMolDuckDB


def run(db_path: str | Path, overwrite: bool = True) -> None:
    db = RetroMolDuckDB.create(db_path, overwrite=overwrite)
    db.close()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    run(db_path=args.db_path, overwrite=args.overwrite)


if __name__ == "__main__":
    main()
