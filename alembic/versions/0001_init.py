"""Version 0001 of the BioNexus database; initial schema."""

from __future__ import annotations

import sqlalchemy as sa
from alembic import op
from pgvector.sqlalchemy import Vector
from sqlalchemy.dialects.postgresql import ARRAY, BIT, JSONB


revision = "0001_init"
down_revision = None
branch_labels = None
depends_on = None


BIOSYNTHETIC_FP_SIZE = 1024


def upgrade() -> None:
    # pgvector extension required for Vector columns and HNSW indexes.
    op.execute("CREATE EXTENSION IF NOT EXISTS vector;")

    # -------------------------------------------------------------------------
    # Shared tables
    # -------------------------------------------------------------------------

    op.create_table(
        "annotation",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column("scheme", sa.String(length=128), nullable=False),
        sa.Column("key", sa.String(length=128), nullable=False),
        sa.Column("value", sa.String(length=512), nullable=False),
        sa.Column("extra", JSONB(), nullable=True),
        sa.UniqueConstraint(
            "scheme",
            "key",
            "value",
            name="ux_annotation_scheme_key_value",
        ),
    )

    op.create_index(
        "ix_annotation_scheme_key_value",
        "annotation",
        ["scheme", "key", "value"],
        unique=False,
    )

    op.create_table(
        "database_reference",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column("database_name", sa.String(length=128), nullable=False),
        sa.Column("database_identifier", sa.String(length=256), nullable=False),
        sa.Column("name", sa.Text(), nullable=True),
        sa.Column("url", sa.Text(), nullable=True),
        sa.Column("extra", JSONB(), nullable=True),
        sa.UniqueConstraint(
            "database_name",
            "database_identifier",
            name="ux_database_reference_database_identifier",
        ),
    )

    op.create_index(
        "ix_database_reference_database_name",
        "database_reference",
        ["database_name"],
        unique=False,
    )
    op.create_index(
        "ix_database_reference_database_identifier",
        "database_reference",
        ["database_identifier"],
        unique=False,
    )

    # -------------------------------------------------------------------------
    # Core entity tables
    # -------------------------------------------------------------------------

    op.create_table(
        "compound",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column("inchikey", sa.String(length=27), nullable=False),
        sa.Column("smiles", sa.Text(), nullable=False),
        sa.Column("canonical_name", sa.Text(), nullable=True),
        sa.Column("formula", sa.String(length=128), nullable=True),
        sa.Column("mol_weight", sa.Float(), nullable=True),
        sa.Column("c_atom_count", sa.Integer(), nullable=True),
        sa.Column("h_atom_count", sa.Integer(), nullable=True),
        sa.Column("n_atom_count", sa.Integer(), nullable=True),
        sa.Column("o_atom_count", sa.Integer(), nullable=True),
        sa.Column("p_atom_count", sa.Integer(), nullable=True),
        sa.Column("s_atom_count", sa.Integer(), nullable=True),
        sa.Column("f_atom_count", sa.Integer(), nullable=True),
        sa.Column("cl_atom_count", sa.Integer(), nullable=True),
        sa.Column("br_atom_count", sa.Integer(), nullable=True),
        sa.Column("i_atom_count", sa.Integer(), nullable=True),
        sa.Column("morgan_fp", BIT(length=2048), nullable=True),
        sa.Column("morgan_fp_radius", sa.Integer(), nullable=True),
        sa.Column("morgan_fp_n_bits", sa.Integer(), nullable=True),
        sa.Column("parser_payload", JSONB(), nullable=True),
        sa.UniqueConstraint("inchikey", name="ux_compound_inchikey"),
    )

    op.create_index(
        "ix_compound_inchikey",
        "compound",
        ["inchikey"],
        unique=False,
    )

    op.create_table(
        "gene_cluster",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column("record_name", sa.Text(), nullable=False),
        sa.Column("file_name", sa.Text(), nullable=False),
        sa.Column("start_bp", sa.BigInteger(), nullable=False),
        sa.Column("end_bp", sa.BigInteger(), nullable=False),
        sa.Column("strand", sa.String(length=8), nullable=True),
        sa.Column("organism", sa.Text(), nullable=True),
        sa.Column("canonical_name", sa.Text(), nullable=True),
        sa.Column("parser_payload", JSONB(), nullable=True),
        sa.CheckConstraint(
            "start_bp <= end_bp",
            name="ck_gene_cluster_start_before_end",
        ),
        sa.UniqueConstraint(
            "record_name",
            "file_name",
            "start_bp",
            "end_bp",
            name="ux_gene_cluster_location",
        ),
    )

    op.create_index(
        "ix_gene_cluster_record_name",
        "gene_cluster",
        ["record_name"],
        unique=False,
    )
    op.create_index(
        "ix_gene_cluster_file_name",
        "gene_cluster",
        ["file_name"],
        unique=False,
    )

    # -------------------------------------------------------------------------
    # Names / aliases
    # -------------------------------------------------------------------------

    op.create_table(
        "compound_name",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column(
            "compound_id",
            sa.BigInteger(),
            sa.ForeignKey("compound.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("name", sa.Text(), nullable=False),
        sa.Column("name_type", sa.String(length=64), nullable=True),
        sa.Column("source", sa.String(length=128), nullable=True),
        sa.UniqueConstraint(
            "compound_id",
            "name",
            "source",
            name="ux_compound_name_compound_name_source",
        ),
    )

    op.create_index(
        "ix_compound_name_compound_id",
        "compound_name",
        ["compound_id"],
        unique=False,
    )
    op.create_index(
        "ix_compound_name_name",
        "compound_name",
        ["name"],
        unique=False,
    )

    op.create_table(
        "gene_cluster_name",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column(
            "gene_cluster_id",
            sa.BigInteger(),
            sa.ForeignKey("gene_cluster.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("name", sa.Text(), nullable=False),
        sa.Column("name_type", sa.String(length=64), nullable=True),
        sa.Column("source", sa.String(length=128), nullable=True),
        sa.UniqueConstraint(
            "gene_cluster_id",
            "name",
            "source",
            name="ux_gene_cluster_name_cluster_name_source",
        ),
    )

    op.create_index(
        "ix_gene_cluster_name_gene_cluster_id",
        "gene_cluster_name",
        ["gene_cluster_id"],
        unique=False,
    )
    op.create_index(
        "ix_gene_cluster_name_name",
        "gene_cluster_name",
        ["name"],
        unique=False,
    )

    # -------------------------------------------------------------------------
    # Retrieval objects
    # -------------------------------------------------------------------------

    op.create_table(
        "retrieval_object",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column(
            "compound_id",
            sa.BigInteger(),
            sa.ForeignKey("compound.id", ondelete="CASCADE"),
            nullable=True,
            unique=True,
        ),
        sa.Column(
            "gene_cluster_id",
            sa.BigInteger(),
            sa.ForeignKey("gene_cluster.id", ondelete="CASCADE"),
            nullable=True,
            unique=True,
        ),
        sa.Column("bag_of_monomers", ARRAY(sa.Text()), nullable=False),
        sa.Column("bag_of_monomers_model", sa.String(length=128), nullable=True),
        sa.Column("coverage", sa.Float(), nullable=True),
        sa.Column("extra", JSONB(), nullable=True),
        sa.CheckConstraint(
            """
            (
                compound_id IS NOT NULL
                AND gene_cluster_id IS NULL
            )
            OR
            (
                compound_id IS NULL
                AND gene_cluster_id IS NOT NULL
            )
            """,
            name="ck_retrieval_object_exactly_one_owner",
        ),
    )

    op.create_index(
        "ix_retrieval_object_compound_id",
        "retrieval_object",
        ["compound_id"],
        unique=False,
    )
    op.create_index(
        "ix_retrieval_object_gene_cluster_id",
        "retrieval_object",
        ["gene_cluster_id"],
        unique=False,
    )
    op.create_index(
        "ix_retrieval_object_bag_of_monomers_gin",
        "retrieval_object",
        ["bag_of_monomers"],
        unique=False,
        postgresql_using="gin",
    )

    # -------------------------------------------------------------------------
    # Primary sequences
    # -------------------------------------------------------------------------

    op.create_table(
        "primary_sequence",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column(
            "retrieval_object_id",
            sa.BigInteger(),
            sa.ForeignKey("retrieval_object.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("source", sa.String(length=128), nullable=False),
        sa.Column("label", sa.String(length=256), nullable=True),
        sa.Column("length", sa.Integer(), nullable=False),
        sa.Column("tokens", ARRAY(sa.Text()), nullable=False),
        sa.Column("biosynthetic_fp", Vector(BIOSYNTHETIC_FP_SIZE), nullable=True),
        sa.Column("biosynthetic_fp_model", sa.String(length=128), nullable=True),
        sa.Column("score", sa.Float(), nullable=True),
        sa.Column("extra", JSONB(), nullable=True),
        sa.CheckConstraint(
            "length >= 0",
            name="ck_primary_sequence_length_nonnegative",
        ),
    )

    op.create_index(
        "ix_primary_sequence_retrieval_object_id",
        "primary_sequence",
        ["retrieval_object_id"],
        unique=False,
    )
    op.create_index(
        "ix_primary_sequence_source",
        "primary_sequence",
        ["source"],
        unique=False,
    )
    op.create_index(
        "ix_primary_sequence_tokens_gin",
        "primary_sequence",
        ["tokens"],
        unique=False,
        postgresql_using="gin",
    )

    op.execute(
        """
        CREATE INDEX IF NOT EXISTS ix_primary_sequence_biosynthetic_fp_hnsw
        ON primary_sequence USING hnsw (biosynthetic_fp vector_cosine_ops);
        """
    )

    op.create_table(
        "primary_sequence_monomer",
        sa.Column("id", sa.BigInteger(), primary_key=True, nullable=False),
        sa.Column(
            "sequence_id",
            sa.BigInteger(),
            sa.ForeignKey("primary_sequence.id", ondelete="CASCADE"),
            nullable=False,
        ),
        sa.Column("position", sa.Integer(), nullable=False),
        sa.Column("token", sa.String(length=128), nullable=False),
        sa.Column("token_class", sa.String(length=128), nullable=True),
        sa.Column("raw_token", sa.String(length=256), nullable=True),
        sa.Column("confidence", sa.Float(), nullable=True),
        sa.Column("extra", JSONB(), nullable=True),
        sa.UniqueConstraint(
            "sequence_id",
            "position",
            name="ux_primary_sequence_monomer_sequence_position",
        ),
        sa.CheckConstraint(
            "position >= 0",
            name="ck_primary_sequence_monomer_position_nonnegative",
        ),
    )

    op.create_index(
        "ix_primary_sequence_monomer_sequence_id",
        "primary_sequence_monomer",
        ["sequence_id"],
        unique=False,
    )
    op.create_index(
        "ix_primary_sequence_monomer_token",
        "primary_sequence_monomer",
        ["token"],
        unique=False,
    )
    op.create_index(
        "ix_primary_sequence_monomer_position_token",
        "primary_sequence_monomer",
        ["position", "token"],
        unique=False,
    )
    op.create_index(
        "ix_primary_sequence_monomer_sequence_position",
        "primary_sequence_monomer",
        ["sequence_id", "position"],
        unique=False,
    )

    # -------------------------------------------------------------------------
    # Association tables
    # -------------------------------------------------------------------------

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
        "gene_cluster_annotation",
        sa.Column(
            "gene_cluster_id",
            sa.BigInteger(),
            sa.ForeignKey("gene_cluster.id", ondelete="CASCADE"),
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
            sa.ForeignKey("database_reference.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
    )

    op.create_table(
        "gene_cluster_reference",
        sa.Column(
            "gene_cluster_id",
            sa.BigInteger(),
            sa.ForeignKey("gene_cluster.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "reference_id",
            sa.BigInteger(),
            sa.ForeignKey("database_reference.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
    )

    op.create_table(
        "annotation_reference",
        sa.Column(
            "annotation_id",
            sa.BigInteger(),
            sa.ForeignKey("annotation.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
        sa.Column(
            "reference_id",
            sa.BigInteger(),
            sa.ForeignKey("database_reference.id", ondelete="CASCADE"),
            primary_key=True,
            nullable=False,
        ),
    )


def downgrade() -> None:
    # Association tables first.
    op.drop_table("annotation_reference")
    op.drop_table("gene_cluster_reference")
    op.drop_table("compound_reference")
    op.drop_table("gene_cluster_annotation")
    op.drop_table("compound_annotation")

    # Primary sequence tables.
    op.drop_index(
        "ix_primary_sequence_monomer_sequence_position",
        table_name="primary_sequence_monomer",
    )
    op.drop_index(
        "ix_primary_sequence_monomer_position_token",
        table_name="primary_sequence_monomer",
    )
    op.drop_index(
        "ix_primary_sequence_monomer_token",
        table_name="primary_sequence_monomer",
    )
    op.drop_index(
        "ix_primary_sequence_monomer_sequence_id",
        table_name="primary_sequence_monomer",
    )
    op.drop_table("primary_sequence_monomer")

    op.execute("DROP INDEX IF EXISTS ix_primary_sequence_biosynthetic_fp_hnsw;")
    op.drop_index("ix_primary_sequence_tokens_gin", table_name="primary_sequence")
    op.drop_index("ix_primary_sequence_source", table_name="primary_sequence")
    op.drop_index("ix_primary_sequence_retrieval_object_id", table_name="primary_sequence")
    op.drop_table("primary_sequence")

    # Retrieval table.
    op.drop_index("ix_retrieval_object_bag_of_monomers_gin", table_name="retrieval_object")
    op.drop_index("ix_retrieval_object_gene_cluster_id", table_name="retrieval_object")
    op.drop_index("ix_retrieval_object_compound_id", table_name="retrieval_object")
    op.drop_table("retrieval_object")

    # Names.
    op.drop_index("ix_gene_cluster_name_name", table_name="gene_cluster_name")
    op.drop_index("ix_gene_cluster_name_gene_cluster_id", table_name="gene_cluster_name")
    op.drop_table("gene_cluster_name")

    op.drop_index("ix_compound_name_name", table_name="compound_name")
    op.drop_index("ix_compound_name_compound_id", table_name="compound_name")
    op.drop_table("compound_name")

    # Core entities.
    op.drop_index("ix_gene_cluster_file_name", table_name="gene_cluster")
    op.drop_index("ix_gene_cluster_record_name", table_name="gene_cluster")
    op.drop_table("gene_cluster")

    op.drop_index("ix_compound_inchikey", table_name="compound")
    op.drop_table("compound")

    # Shared tables.
    op.drop_index(
        "ix_database_reference_database_identifier",
        table_name="database_reference",
    )
    op.drop_index(
        "ix_database_reference_database_name",
        table_name="database_reference",
    )
    op.drop_table("database_reference")

    op.drop_index("ix_annotation_scheme_key_value", table_name="annotation")
    op.drop_table("annotation")