"""Database engine and session setup for BioNexus."""

import os

from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker


DB_URL = os.getenv("BIONEXUS_DB_URL")
if DB_URL is None:
    raise ValueError("BIONEXUS_DB_URL environment variable is not set")
engine = create_engine(DB_URL, pool_pre_ping=True)
SessionLocal = sessionmaker(bind=engine, autoflush=False, expire_on_commit=False)