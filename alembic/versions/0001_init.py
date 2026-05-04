"""Version 0001 of the RetroMol database; initial schema."""

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import JSONB
from pgvector.sqlalchemy import Vector
from alembic import op


revision = "0001_init"
down_revision = None
branch_labels = None
depends_on = None


def upgrade():
    # pgvector extension (required for Vector columns + ANN indexes)
    op.execute("CREATE EXTENSION IF NOT EXISTS vector;")

    op.create_table(
        "annotation",
        sa.Column("id", sa.BigInteger, primary_key=True, nullable=False),
        sa.Column("scheme", sa.String(64), nullable=False),
        sa.Column("key", sa.String(64), nullable=False),
        sa.Column("value", sa.String(256), nullable=False),
        sa.UniqueConstraint(
            "scheme",
            "key",
            "value",
            name="ux_annotation_scheme_key_value",
        ),
    )

    op.create_table(
        "compound",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column("inchikey", sa.String(length=27), nullable=False, unique=True),
        sa.Column("smiles", sa.Text(), nullable=False),
        sa.Column("mol_weight", sa.Float(), nullable=False),
        sa.Column("c_atom_count", sa.Integer(), nullable=False),
        sa.Column("h_atom_count", sa.Integer(), nullable=False),
        sa.Column("n_atom_count", sa.Integer(), nullable=False),
        sa.Column("o_atom_count", sa.Integer(), nullable=False),
        sa.Column("p_atom_count", sa.Integer(), nullable=False),
        sa.Column("s_atom_count", sa.Integer(), nullable=False),
        sa.Column("f_atom_count", sa.Integer(), nullable=False),
        sa.Column("cl_atom_count", sa.Integer(), nullable=False),
        sa.Column("br_atom_count", sa.Integer(), nullable=False),
        sa.Column("i_atom_count", sa.Integer(), nullable=False),
        sa.Column("morgan_fp", Vector(2048), nullable=False),
        sa.Column("retromol_fp_counted", Vector(1024), nullable=False),
        sa.Column("retromol_fp_binary", Vector(1024), nullable=False),
        sa.Column("retromol", JSONB(), nullable=False),
        sa.Column("coverage", sa.Float(), nullable=False),
    )

    op.create_index("ix_compound_inchikey", "compound", ["inchikey"], unique=False)

    op.create_table(
        "candidate_cluster",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column("record_name", sa.Text(), nullable=False),
        sa.Column("file_name", sa.Text(), nullable=False),
        sa.Column("start_bp", sa.BigInteger(), nullable=False),
        sa.Column("end_bp", sa.BigInteger(), nullable=False),
        sa.Column("retromol_fp_counted_by_orf", Vector(1024), nullable=False),
        sa.Column("retromol_fp_binary_by_orf", Vector(1024), nullable=False),
        sa.Column("retromol_fp_counted_by_region", Vector(1024), nullable=False),
        sa.Column("retromol_fp_binary_by_region", Vector(1024), nullable=False),
        sa.Column("biocracker", JSONB(), nullable=False),
        sa.UniqueConstraint(
            "record_name",
            "file_name",
            "start_bp",
            "end_bp",
            name="ux_candidate_cluster_recordname_filename_start_end",
        ),
    )

    op.create_table(
        "reference",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column("name", sa.String(length=256), nullable=False),
        sa.Column("database_name", sa.String(length=128), nullable=False),
        sa.Column("database_identifier", sa.String(length=128), nullable=False),
        sa.UniqueConstraint(
            "name",
            "database_name",
            "database_identifier",
            name="ux_reference_name_dbname_dbid",
        ),
    )

    op.create_table(
        "compound_annotation",
        sa.Column(
            "compound_id",
            sa.BigInteger(),
            sa.ForeignKey("compound.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "annotation_id",
            sa.BigInteger(),
            sa.ForeignKey("annotation.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
    )

    op.create_table(
        "candidate_cluster_annotation",
        sa.Column(
            "candidate_cluster_id",
            sa.BigInteger(),
            sa.ForeignKey("candidate_cluster.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "annotation_id",
            sa.BigInteger(),
            sa.ForeignKey("annotation.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
    )

    op.create_table(
        "compound_reference",
        sa.Column(
            "compound_id",
            sa.BigInteger(),
            sa.ForeignKey("compound.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "reference_id",
            sa.BigInteger(),
            sa.ForeignKey("reference.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
    )

    op.create_table(
        "candidate_cluster_reference",
        sa.Column(
            "candidate_cluster_id",
            sa.BigInteger(),
            sa.ForeignKey("candidate_cluster.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "reference_id",
            sa.BigInteger(),
            sa.ForeignKey("reference.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
    )

    # pgvector ANN indexes (cosine)
    op.execute("""
        CREATE INDEX IF NOT EXISTS ix_compound_retromol_fp_counted_hnsw
        ON compound USING hnsw (retromol_fp_counted vector_cosine_ops);
    """)
    op.execute("""
        CREATE INDEX IF NOT EXISTS ix_compound_retromol_fp_binary_hnsw
        ON compound USING hnsw (retromol_fp_binary vector_cosine_ops);
    """)
    op.execute("""
        CREATE INDEX IF NOT EXISTS ix_candidate_cluster_retromol_fp_counted_by_orf_hnsw
        ON candidate_cluster USING hnsw (retromol_fp_counted_by_orf vector_cosine_ops);
    """)
    op.execute("""
        CREATE INDEX IF NOT EXISTS ix_candidate_cluster_retromol_fp_binary_by_orf_hnsw
        ON candidate_cluster USING hnsw (retromol_fp_binary_by_orf vector_cosine_ops);
    """)
    op.execute("""
        CREATE INDEX IF NOT EXISTS ix_candidate_cluster_retromol_fp_counted_by_region_hnsw
        ON candidate_cluster USING hnsw (retromol_fp_counted_by_region vector_cosine_ops);
    """)
    op.execute("""
        CREATE INDEX IF NOT EXISTS ix_candidate_cluster_retromol_fp_binary_by_region_hnsw
        ON candidate_cluster USING hnsw (retromol_fp_binary_by_region vector_cosine_ops);
    """)


def downgrade():
    # Drop ANN indexes
    op.execute("DROP INDEX IF EXISTS ix_candidate_cluster_retromol_fp_binary_by_region_hnsw;")
    op.execute("DROP INDEX IF EXISTS ix_candidate_cluster_retromol_fp_counted_by_region_hnsw;")
    op.execute("DROP INDEX IF EXISTS ix_candidate_cluster_retromol_fp_binary_by_orf_hnsw;")
    op.execute("DROP INDEX IF EXISTS ix_candidate_cluster_retromol_fp_counted_by_orf_hnsw;")
    op.execute("DROP INDEX IF EXISTS ix_compound_retromol_fp_binary_hnsw;")
    op.execute("DROP INDEX IF EXISTS ix_compound_retromol_fp_counted_hnsw;")

    # Drop tables (association tables first due to FKs)
    op.drop_table("candidate_cluster_reference")
    op.drop_table("compound_reference")

    op.drop_table("candidate_cluster_annotation")
    op.drop_table("compound_annotation")

    op.drop_table("reference")

    op.drop_table("candidate_cluster")

    op.drop_index("ix_compound_inchikey", table_name="compound")
    op.drop_table("compound")

    op.drop_table("annotation")