import os
from pathlib import Path


from retromol_database.duckdb import RetroMolDuckDB, open_database


def duckdb_path_from_env() -> Path:
    path = os.getenv("RETROMOL_DUCKDB_PATH")
    if not path:
        raise Exception("RETROMOL_DUCKDB_PATH is not set")
    return Path(path).expanduser()


def open_retromol_db() -> RetroMolDuckDB:
    return open_database(duckdb_path_from_env(), read_only=True)


def check_database_ready() -> None:
    path = duckdb_path_from_env()
    if not path.exists():
        raise FileNotFoundError(path)

    with open_retromol_db() as db:
        db.count()
