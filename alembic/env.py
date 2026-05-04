"""Alembic migration environment configuration file."""

import os
from logging.config import fileConfig
from pathlib import Path

from sqlalchemy import engine_from_config, pool
from alembic import context

from retromol_database.db import models as bnx_models

# Make sure we import local repo version, not installed package
REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT) not in os.sys.path:
    os.sys.path.insert(0, str(REPO_ROOT))

# Load .env from repo root
try:
    from dotenv import load_dotenv
    load_dotenv(Path(__file__).resolve().parents[1] / ".env")
except Exception:
    pass

config = context.config
if config.config_file_name is not None:
    fileConfig(config.config_file_name)

# Import module that defines base and all models
target_metadata = bnx_models.Base.metadata


def get_url() -> str:
    """
    Get the database URL from environment variable.

    :return: database URL string
    :raises RuntimeError: if BIONEXUS_DB_URL is not set
    """
    url = os.getenv("BIONEXUS_DB_URL")
    if not url:
        raise RuntimeError("BIONEXUS_DB_URL environment variable not set")
    return url


def run_migrations_offline() -> None:
    """
    Run migrations in 'offline' mode.
    """
    url = get_url()
    context.configure(
        url=url,
        target_metadata=target_metadata,
        literal_binds=True,
        compare_type=True,
        compare_server_default=True,
    )
    with context.begin_transaction():
        context.run_migrations()


def run_migrations_online() -> None:
    """
    Run migrations in 'online' mode.
    """
    connectable = engine_from_config({"sqlalchemy.url": get_url()}, prefix="sqlalchemy.", poolclass=pool.NullPool)
    with connectable.connect() as connection:
        context.configure(
            connection=connection,
            target_metadata=target_metadata,
            compare_type=True,
            compare_server_default=True,
        )
        with context.begin_transaction():
            context.run_migrations()


# Determine mode and run migrations accordingly
if context.is_offline_mode():
    run_migrations_offline()
else:
    run_migrations_online()