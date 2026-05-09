"""Database models for BioNexus."""

from __future__ import annotations

from typing import Any

import sqlalchemy as sa
from pgvector.sqlalchemy import Vector
from sqlalchemy.dialects.postgresql import ARRAY, BIT, JSONB
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship


BIOSYNTHETIC_FP_SIZE = 1024


class Base(DeclarativeBase):
    """Base class for all database models."""

    pass


# =============================================================================
# Association tables
# =============================================================================


compound_annotation = sa.Table(
    "compound_annotation",
    Base.metadata,
    sa.Column(
        "compound_id",
        sa.BigInteger,
        sa.ForeignKey("compound.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    sa.Column(
        "annotation_id",
        sa.BigInteger,
        sa.ForeignKey("annotation.id", ondelete="CASCADE"),
        primary_key=True,
    ),
)


gene_cluster_annotation = sa.Table(
    "gene_cluster_annotation",
    Base.metadata,
    sa.Column(
        "gene_cluster_id",
        sa.BigInteger,
        sa.ForeignKey("gene_cluster.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    sa.Column(
        "annotation_id",
        sa.BigInteger,
        sa.ForeignKey("annotation.id", ondelete="CASCADE"),
        primary_key=True,
    ),
)


compound_reference = sa.Table(
    "compound_reference",
    Base.metadata,
    sa.Column(
        "compound_id",
        sa.BigInteger,
        sa.ForeignKey("compound.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    sa.Column(
        "reference_id",
        sa.BigInteger,
        sa.ForeignKey("database_reference.id", ondelete="CASCADE"),
        primary_key=True,
    ),
)


gene_cluster_reference = sa.Table(
    "gene_cluster_reference",
    Base.metadata,
    sa.Column(
        "gene_cluster_id",
        sa.BigInteger,
        sa.ForeignKey("gene_cluster.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    sa.Column(
        "reference_id",
        sa.BigInteger,
        sa.ForeignKey("database_reference.id", ondelete="CASCADE"),
        primary_key=True,
    ),
)


annotation_reference = sa.Table(
    "annotation_reference",
    Base.metadata,
    sa.Column(
        "annotation_id",
        sa.BigInteger,
        sa.ForeignKey("annotation.id", ondelete="CASCADE"),
        primary_key=True,
    ),
    sa.Column(
        "reference_id",
        sa.BigInteger,
        sa.ForeignKey("database_reference.id", ondelete="CASCADE"),
        primary_key=True,
    ),
)


# =============================================================================
# Shared models
# =============================================================================


class DatabaseReference(Base):
    """Database or literature reference.

    Can be linked to compounds, gene clusters, and annotations.
    """

    __tablename__ = "database_reference"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    database_name: Mapped[str] = mapped_column(sa.String(128), nullable=False)
    database_identifier: Mapped[str] = mapped_column(sa.String(256), nullable=False)

    name: Mapped[str | None] = mapped_column(sa.Text, nullable=True)
    url: Mapped[str | None] = mapped_column(sa.Text, nullable=True)
    extra: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    __table_args__ = (
        sa.UniqueConstraint(
            "database_name",
            "database_identifier",
            name="ux_database_reference_database_identifier",
        ),
        sa.Index("ix_database_reference_database_name", "database_name"),
        sa.Index("ix_database_reference_database_identifier", "database_identifier"),
    )

    compounds: Mapped[list[Compound]] = relationship(
        "Compound",
        secondary=compound_reference,
        back_populates="references",
        lazy="selectin",
    )

    gene_clusters: Mapped[list[GeneCluster]] = relationship(
        "GeneCluster",
        secondary=gene_cluster_reference,
        back_populates="references",
        lazy="selectin",
    )

    annotations: Mapped[list[Annotation]] = relationship(
        "Annotation",
        secondary=annotation_reference,
        back_populates="references",
        lazy="selectin",
    )


class Annotation(Base):
    """Annotation used for enrichment."""

    __tablename__ = "annotation"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    scheme: Mapped[str] = mapped_column(sa.String(128), nullable=False)
    key: Mapped[str] = mapped_column(sa.String(128), nullable=False)
    value: Mapped[str] = mapped_column(sa.String(512), nullable=False)

    extra: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    __table_args__ = (
        sa.UniqueConstraint(
            "scheme",
            "key",
            "value",
            name="ux_annotation_scheme_key_value",
        ),
        sa.Index("ix_annotation_scheme_key_value", "scheme", "key", "value"),
    )

    compounds: Mapped[list[Compound]] = relationship(
        "Compound",
        secondary=compound_annotation,
        back_populates="annotations",
        lazy="selectin",
    )

    gene_clusters: Mapped[list[GeneCluster]] = relationship(
        "GeneCluster",
        secondary=gene_cluster_annotation,
        back_populates="annotations",
        lazy="selectin",
    )

    references: Mapped[list[DatabaseReference]] = relationship(
        "DatabaseReference",
        secondary=annotation_reference,
        back_populates="annotations",
        lazy="selectin",
    )


# =============================================================================
# Core entity models
# =============================================================================


class Compound(Base):
    """Chemical compound.

    A compound can have:
    - zero-to-many names
    - zero-to-many database references
    - zero-to-many annotations
    - zero-or-one retrieval object
    """

    __tablename__ = "compound"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    inchikey: Mapped[str] = mapped_column(sa.String(27), nullable=False, unique=True, index=True)
    smiles: Mapped[str] = mapped_column(sa.Text, nullable=False)

    canonical_name: Mapped[str | None] = mapped_column(sa.Text, nullable=True)
    formula: Mapped[str | None] = mapped_column(sa.String(128), nullable=True)
    mol_weight: Mapped[float | None] = mapped_column(sa.Float, nullable=True)

    c_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    h_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    n_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    o_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    p_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    s_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    f_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    cl_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    br_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    i_atom_count: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)

    # Normal Morgan fingerprint, not pgvector, not HNSW.
    # Store RDKit BitVect as a 2048-character bit string.
    morgan_fp: Mapped[str | None] = mapped_column(BIT(length=2048), nullable=True)
    morgan_fp_radius: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)
    morgan_fp_n_bits: Mapped[int | None] = mapped_column(sa.Integer, nullable=True)

    # Raw parser output, e.g. RetroMol output.
    parser_payload: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    names: Mapped[list[CompoundName]] = relationship(
        "CompoundName",
        back_populates="compound",
        cascade="all, delete-orphan",
        lazy="selectin",
    )

    references: Mapped[list[DatabaseReference]] = relationship(
        "DatabaseReference",
        secondary=compound_reference,
        back_populates="compounds",
        lazy="selectin",
    )

    annotations: Mapped[list[Annotation]] = relationship(
        "Annotation",
        secondary=compound_annotation,
        back_populates="compounds",
        lazy="selectin",
    )

    retrieval: Mapped[RetrievalObject | None] = relationship(
        "RetrievalObject",
        back_populates="compound",
        uselist=False,
        cascade="all, delete-orphan",
        lazy="selectin",
    )


class GeneCluster(Base):
    """Biosynthetic gene cluster."""

    __tablename__ = "gene_cluster"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    record_name: Mapped[str] = mapped_column(sa.Text, nullable=False)
    file_name: Mapped[str] = mapped_column(sa.Text, nullable=False)
    start_bp: Mapped[int] = mapped_column(sa.BigInteger, nullable=False)
    end_bp: Mapped[int] = mapped_column(sa.BigInteger, nullable=False)

    strand: Mapped[str | None] = mapped_column(sa.String(8), nullable=True)
    organism: Mapped[str | None] = mapped_column(sa.Text, nullable=True)
    canonical_name: Mapped[str | None] = mapped_column(sa.Text, nullable=True)

    # Raw parser output, e.g. antiSMASH/BioCracker output.
    parser_payload: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    __table_args__ = (
        sa.UniqueConstraint(
            "record_name",
            "file_name",
            "start_bp",
            "end_bp",
            name="ux_gene_cluster_location",
        ),
        sa.Index("ix_gene_cluster_record_name", "record_name"),
        sa.Index("ix_gene_cluster_file_name", "file_name"),
        sa.CheckConstraint("start_bp <= end_bp", name="ck_gene_cluster_start_before_end"),
    )

    names: Mapped[list[GeneClusterName]] = relationship(
        "GeneClusterName",
        back_populates="gene_cluster",
        cascade="all, delete-orphan",
        lazy="selectin",
    )

    references: Mapped[list[DatabaseReference]] = relationship(
        "DatabaseReference",
        secondary=gene_cluster_reference,
        back_populates="gene_clusters",
        lazy="selectin",
    )

    annotations: Mapped[list[Annotation]] = relationship(
        "Annotation",
        secondary=gene_cluster_annotation,
        back_populates="gene_clusters",
        lazy="selectin",
    )

    retrieval: Mapped[RetrievalObject | None] = relationship(
        "RetrievalObject",
        back_populates="gene_cluster",
        uselist=False,
        cascade="all, delete-orphan",
        lazy="selectin",
    )


# =============================================================================
# Names / aliases
# =============================================================================


class CompoundName(Base):
    """Alternative compound name."""

    __tablename__ = "compound_name"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    compound_id: Mapped[int] = mapped_column(
        sa.BigInteger,
        sa.ForeignKey("compound.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )

    name: Mapped[str] = mapped_column(sa.Text, nullable=False)
    name_type: Mapped[str | None] = mapped_column(sa.String(64), nullable=True)
    source: Mapped[str | None] = mapped_column(sa.String(128), nullable=True)

    __table_args__ = (
        sa.UniqueConstraint(
            "compound_id",
            "name",
            "source",
            name="ux_compound_name_compound_name_source",
        ),
        sa.Index("ix_compound_name_name", "name"),
    )

    compound: Mapped[Compound] = relationship(
        "Compound",
        back_populates="names",
    )


class GeneClusterName(Base):
    """Alternative BGC/gene-cluster name."""

    __tablename__ = "gene_cluster_name"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    gene_cluster_id: Mapped[int] = mapped_column(
        sa.BigInteger,
        sa.ForeignKey("gene_cluster.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )

    name: Mapped[str] = mapped_column(sa.Text, nullable=False)
    name_type: Mapped[str | None] = mapped_column(sa.String(64), nullable=True)
    source: Mapped[str | None] = mapped_column(sa.String(128), nullable=True)

    __table_args__ = (
        sa.UniqueConstraint(
            "gene_cluster_id",
            "name",
            "source",
            name="ux_gene_cluster_name_cluster_name_source",
        ),
        sa.Index("ix_gene_cluster_name_name", "name"),
    )

    gene_cluster: Mapped[GeneCluster] = relationship(
        "GeneCluster",
        back_populates="names",
    )


# =============================================================================
# Retrieval object
# =============================================================================


class RetrievalObject(Base):
    """Retrieval profile for either a compound or a gene cluster.

    A retrieval object owns:
    - an unordered bag of identified monomers
    - zero-to-many primary sequences
    """

    __tablename__ = "retrieval_object"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    compound_id: Mapped[int | None] = mapped_column(
        sa.BigInteger,
        sa.ForeignKey("compound.id", ondelete="CASCADE"),
        nullable=True,
        unique=True,
    )

    gene_cluster_id: Mapped[int | None] = mapped_column(
        sa.BigInteger,
        sa.ForeignKey("gene_cluster.id", ondelete="CASCADE"),
        nullable=True,
        unique=True,
    )

    # Unordered bag of all identified monomers for this compound or gene cluster.
    # Keep duplicates if you want counted-bag semantics, e.g. ["Leu", "Leu", "Val"].
    # Use array containment operators for simple membership queries.
    bag_of_monomers: Mapped[list[str]] = mapped_column(
        ARRAY(sa.Text),
        nullable=False,
        default=list,
    )

    bag_of_monomers_model: Mapped[str | None] = mapped_column(sa.String(128), nullable=True)

    coverage: Mapped[float | None] = mapped_column(sa.Float, nullable=True)

    extra: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    __table_args__ = (
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
        sa.Index("ix_retrieval_object_compound_id", "compound_id"),
        sa.Index("ix_retrieval_object_gene_cluster_id", "gene_cluster_id"),
        sa.Index("ix_retrieval_object_bag_of_monomers_gin", "bag_of_monomers", postgresql_using="gin"),
    )

    compound: Mapped[Compound | None] = relationship(
        "Compound",
        back_populates="retrieval",
    )

    gene_cluster: Mapped[GeneCluster | None] = relationship(
        "GeneCluster",
        back_populates="retrieval",
    )

    primary_sequences: Mapped[list[PrimarySequence]] = relationship(
        "PrimarySequence",
        back_populates="retrieval",
        cascade="all, delete-orphan",
        lazy="selectin",
    )


# =============================================================================
# Primary sequences
# =============================================================================


class PrimarySequence(Base):
    """Linear sequence of building blocks.

    Belongs to one retrieval object.

    Each primary sequence has its own fingerprint so you can first retrieve
    candidates through sequence-level vector similarity and then refine with
    explicit monomer-position queries.
    """

    __tablename__ = "primary_sequence"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    retrieval_object_id: Mapped[int] = mapped_column(
        sa.BigInteger,
        sa.ForeignKey("retrieval_object.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )

    source: Mapped[str] = mapped_column(sa.String(128), nullable=False)
    label: Mapped[str | None] = mapped_column(sa.String(256), nullable=True)

    length: Mapped[int] = mapped_column(sa.Integer, nullable=False)

    # Full sequence for convenience.
    tokens: Mapped[list[str]] = mapped_column(ARRAY(sa.Text), nullable=False)

    # Sequence-specific fingerprint.
    biosynthetic_fp: Mapped[list[float] | None] = mapped_column(
        Vector(BIOSYNTHETIC_FP_SIZE),
        nullable=True,
    )
    biosynthetic_fp_model: Mapped[str | None] = mapped_column(sa.String(128), nullable=True)

    score: Mapped[float | None] = mapped_column(sa.Float, nullable=True)

    extra: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    __table_args__ = (
        sa.CheckConstraint("length >= 0", name="ck_primary_sequence_length_nonnegative"),
        sa.Index("ix_primary_sequence_tokens_gin", "tokens", postgresql_using="gin"),
        sa.Index("ix_primary_sequence_source", "source"),
        sa.Index(
            "ix_primary_sequence_biosynthetic_fp_hnsw",
            "biosynthetic_fp",
            postgresql_using="hnsw",
            postgresql_ops={"biosynthetic_fp": "vector_cosine_ops"},
        ),
    )

    retrieval: Mapped[RetrievalObject] = relationship(
        "RetrievalObject",
        back_populates="primary_sequences",
    )

    monomers: Mapped[list[PrimarySequenceMonomer]] = relationship(
        "PrimarySequenceMonomer",
        back_populates="sequence",
        cascade="all, delete-orphan",
        order_by="PrimarySequenceMonomer.position",
        lazy="selectin",
    )


class PrimarySequenceMonomer(Base):
    """One monomer/building block at one position in a primary sequence."""

    __tablename__ = "primary_sequence_monomer"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    sequence_id: Mapped[int] = mapped_column(
        sa.BigInteger,
        sa.ForeignKey("primary_sequence.id", ondelete="CASCADE"),
        nullable=False,
        index=True,
    )

    # Zero-based position.
    position: Mapped[int] = mapped_column(sa.Integer, nullable=False)

    # Normalized searchable token, e.g. "Leu", "A1", "B2".
    token: Mapped[str] = mapped_column(sa.String(128), nullable=False)

    token_class: Mapped[str | None] = mapped_column(sa.String(128), nullable=True)
    raw_token: Mapped[str | None] = mapped_column(sa.String(256), nullable=True)
    confidence: Mapped[float | None] = mapped_column(sa.Float, nullable=True)

    extra: Mapped[dict[str, Any] | None] = mapped_column(JSONB, nullable=True)

    __table_args__ = (
        sa.UniqueConstraint(
            "sequence_id",
            "position",
            name="ux_primary_sequence_monomer_sequence_position",
        ),
        sa.CheckConstraint("position >= 0", name="ck_primary_sequence_monomer_position_nonnegative"),
        sa.Index("ix_primary_sequence_monomer_token", "token"),
        sa.Index("ix_primary_sequence_monomer_position_token", "position", "token"),
        sa.Index("ix_primary_sequence_monomer_sequence_position", "sequence_id", "position"),
    )

    sequence: Mapped[PrimarySequence] = relationship(
        "PrimarySequence",
        back_populates="monomers",
    )