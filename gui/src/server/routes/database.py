"""Database connection setup using SQLAlchemy."""

import os

from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine


def dsn_from_env() -> str:
    """
    Construct the Postgres DSN from environment variables.

    :return: the Postgres DSN string
    """
    dsn = os.getenv("DATABASE_URL")
    if dsn:
        # If plain postgresql:// URL was provided, force psycopg v3 driver
        if dsn.startswith("postgresql://"):
            dsn = dsn.replace("postgresql://", "postgresql+psycopg://", 1)
        # When explicitly using psycopg2, upgrade it too (optional but helpful)
        if dsn.startswith("postgresql+psycopg2://"):
            dsn = dsn.replace("postgresql+psycopg2://", "postgresql+psycopg://", 1)
        return dsn

    host = os.getenv("DB_HOST", "db")
    port = os.getenv("DB_PORT", "5432")
    name = os.getenv("DB_NAME", "bionexus")
    user = os.getenv("DB_USER", "app_ro")
    pwd = os.getenv("DB_PASS") or os.getenv("DB_PASSWORD", "apppass_ro")

    url = f"postgresql+psycopg://{user}:{pwd}@{host}:{port}/{name}"
    print(f"Constructed DSN: {url}")
    return url


engine = create_engine(dsn_from_env(), pool_pre_ping=True)
SessionLocal = sessionmaker(bind=engine, autoflush=False, expire_on_commit=False)