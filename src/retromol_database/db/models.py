"""Database models for BioNexus."""

import sqlalchemy as sa
from sqlalchemy.dialects.postgresql import ARRAY, JSONB
from sqlalchemy.orm import DeclarativeBase, Mapped, mapped_column, relationship
from pgvector.sqlalchemy import Vector


class Base(DeclarativeBase):
    """
    Base class for all database models.
    """

    pass


compound_annotation = sa.Table(
    "compound_annotation",
    Base.metadata,
    sa.Column("compound_id", sa.BigInteger, sa.ForeignKey("compound.id", ondelete="CASCADE"), primary_key=True),
    sa.Column("annotation_id", sa.BigInteger, sa.ForeignKey("annotation.id", ondelete="CASCADE"), primary_key=True),
)

candidate_cluster_annotation = sa.Table(
    "candidate_cluster_annotation",
    Base.metadata,
    sa.Column("candidate_cluster_id", sa.BigInteger, sa.ForeignKey("candidate_cluster.id", ondelete="CASCADE"),
              primary_key=True),
    sa.Column("annotation_id", sa.BigInteger, sa.ForeignKey("annotation.id", ondelete="CASCADE"), primary_key=True),
)

compound_reference = sa.Table(
    "compound_reference",
    Base.metadata,
    sa.Column("compound_id", sa.BigInteger, sa.ForeignKey("compound.id", ondelete="CASCADE"), primary_key=True),
    sa.Column("reference_id", sa.BigInteger, sa.ForeignKey("reference.id", ondelete="CASCADE"), primary_key=True),
)

candidate_cluster_reference = sa.Table(
    "candidate_cluster_reference",
    Base.metadata,
    sa.Column("candidate_cluster_id", sa.BigInteger, sa.ForeignKey("candidate_cluster.id", ondelete="CASCADE"),
              primary_key=True),
    sa.Column("reference_id", sa.BigInteger, sa.ForeignKey("reference.id", ondelete="CASCADE"), primary_key=True),
)


class Reference(Base):
    """
    Describes a literature or database reference linked to either or both a Compound and/or a CandidateCluster.
    """

    __tablename__ = "reference"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    name: Mapped[str] = mapped_column(sa.String(256), nullable=False)
    database_name: Mapped[str] = mapped_column(sa.String(128), nullable=False)
    database_identifier: Mapped[str] = mapped_column(sa.String(128), nullable=False)

    # Only one reference can exist for a database per identifier (e.g., only one NPA000001 in NPAtlas)
    __table_args__ = (
        sa.UniqueConstraint(
            "name",
            "database_name",
            "database_identifier",
            name="ux_reference_name_dbname_dbid",
        ),
    )

    compounds = relationship(
        "Compound",
        secondary=compound_reference,
        back_populates="references",
        lazy="selectin",
    )
    candidate_clusters = relationship(
        "CandidateCluster",
        secondary=candidate_cluster_reference,
        back_populates="references",
        lazy="selectin",
    )


class Annotation(Base):
    """
    Describes an annotation linked to either or both a Compound and/or a CandidateCluster.
    """

    __tablename__ = "annotation"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    # Label triplet
    scheme: Mapped[str] = mapped_column(sa.String(64), nullable=False)
    key: Mapped[str] = mapped_column(sa.String(64), nullable=False)
    value: Mapped[str] = mapped_column(sa.String(256), nullable=False)

    __table_args__ = (
        sa.UniqueConstraint(
            "scheme",
            "key",
            "value",
            name="ux_annotation_scheme_key_value",
        ),
    )

    compounds = relationship(
        "Compound",
        secondary=compound_annotation,
        back_populates="annotations",
        lazy="selectin",
    )
    candidate_clusters = relationship(
        "CandidateCluster",
        secondary=candidate_cluster_annotation,
        back_populates="annotations",
        lazy="selectin",
    )


class Compound(Base):
    """
    Describes a chemical structure.
    """

    __tablename__ = "compound"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    # Canonical identifiers/structure
    inchikey: Mapped[str] = mapped_column(sa.String(27), nullable=False, unique=True, index=True)
    smiles: Mapped[str] = mapped_column(sa.Text, nullable=False)

    # Canonical properties
    mol_weight: Mapped[float] = mapped_column(sa.Float, nullable=False)

    # Elemental counts
    c_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    h_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    n_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    o_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    p_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    s_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    f_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    cl_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    br_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)
    i_atom_count: Mapped[int] = mapped_column(sa.Integer, nullable=False)

    # Fingerprints
    morgan_fp: Mapped[list[float] | None] = mapped_column(Vector(2048), nullable=False)
    retromol_fp_binary: Mapped[list[float] | None] = mapped_column(Vector(1024), nullable=False)
    retromol_fp_counted: Mapped[list[float] | None] = mapped_column(Vector(1024), nullable=False)

    # RetroMol parsing results (as serialized dict/JSON)
    retromol: Mapped[dict[str, str] | None] = mapped_column(JSONB, nullable=False)
    coverage: Mapped[float] = mapped_column(sa.Float, nullable=False)

    annotations = relationship(
        "Annotation",
        secondary=compound_annotation,
        back_populates="compounds",
        lazy="selectin",
    )
    references = relationship(
        "Reference",
        secondary=compound_reference,
        back_populates="compounds",
        lazy="selectin",
    )


class CandidateCluster(Base):
    """
    Describes a candidate BGC cluster.
    """

    __tablename__ = "candidate_cluster"

    id: Mapped[int] = mapped_column(sa.BigInteger, primary_key=True)

    # Candidate cluster identifiers
    record_name: Mapped[str] = mapped_column(sa.Text, nullable=False)
    file_name: Mapped[str] = mapped_column(sa.Text, nullable=False)
    start_bp: Mapped[int] = mapped_column(sa.BigInteger, nullable=False)
    end_bp: Mapped[int] = mapped_column(sa.BigInteger, nullable=False)

    # Fingerprints
    retromol_fp_binary_by_orf: Mapped[list[float] | None] = mapped_column(Vector(1024), nullable=False)
    retromol_fp_counted_by_orf: Mapped[list[float] | None] = mapped_column(Vector(1024), nullable=False)
    retromol_fp_binary_by_region: Mapped[list[float] | None] = mapped_column(Vector(1024), nullable=False)
    retromol_fp_counted_by_region: Mapped[list[float] | None] = mapped_column(Vector(1024), nullable=False)

    # BioCracker parsing results (as serialized dict/JSON)
    biocracker: Mapped[dict[str, str] | None] = mapped_column(JSONB, nullable=False)

    __table_args__ = (
        sa.UniqueConstraint(
            "record_name",
            "file_name",
            "start_bp",
            "end_bp",
            name="ux_candidate_cluster_recordname_filename_start_end",
        ),
    )

    annotations = relationship(
        "Annotation",
        secondary=candidate_cluster_annotation,
        back_populates="candidate_clusters",
        lazy="selectin",
    )
    references = relationship(
        "Reference",
        secondary=candidate_cluster_reference,
        back_populates="candidate_clusters",
        lazy="selectin",
    )