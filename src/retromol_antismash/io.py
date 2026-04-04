"""Module containing functionalities for reading and parsing antiSMASH GeneBank files."""

import logging
from pathlib import Path
from dataclasses import dataclass, field
from typing import Literal, Union, TypeAlias, Any, Iterable, TYPE_CHECKING

try:
    from Bio import SeqIO
    BIOPYTHON_AVAILABLE = True
except ImportError:
    SeqIO = None
    BIOPYTHON_AVAILABLE = False

if TYPE_CHECKING:
    from Bio.SeqFeature import FeatureLocation, SeqFeature
    from Bio.SeqRecord import SeqRecord

from retromol_antismash.model import Region, Gene, Strand, Domain


log = logging.getLogger(__name__)


AntismashReadoutLevel: TypeAlias = Literal["region", "cand_cluster"]


@dataclass(frozen=True)
class AntiSmashOptions:
    """
    Options for loading biosynthetic regions from antiSMASH GenBank files.
    
    :param source: Source type.
    :param readout_level: Parsing level for biosynthetic regions.
    :param wanted_qualifiers: List of feature qualifiers to extract from regions.
    :param gene_identifiers: List of gene feature types to look for.
    :param domain_identifiers: List of domain feature types to look for.
    """

    source: Literal["antismash_gbk"] = "antismash_gbk"
    readout_level: AntismashReadoutLevel = "region"
    wanted_qualifiers: list[str] | None = field(default_factory=lambda: ["product"])
    gene_identifiers: list[str] | None = field(default_factory=lambda: ["CDS"])
    domain_identifiers: list[str] | None = field(default_factory=lambda: ["aSDomain"])

    def __post_init__(self) -> None:
        """
        Validate options after initialization.
        """
        valid_levels = {"region", "cand_cluster"}
        if self.readout_level not in valid_levels:
            raise ValueError(f"expected readout_level to be one of {valid_levels}, got {self.readout_level}")


RegionLoadOptions = Union[AntiSmashOptions]


def _iter_regions(record: "SeqRecord", readout_level: AntismashReadoutLevel) -> list["SeqFeature"]:
    """
    Iterate over region features in a Biopython SeqRecord.
    
    :param record: Biopython SeqRecord object.
    :return: List of region SeqFeature objects.
    """
    return [f for f in (record.features or []) if f.type == readout_level]


def _iter_cds(record: "SeqRecord", gene_identifiers: list[str] | None = None) -> list["SeqFeature"]:
    """
    Iterate over CDS features in a Biopython SeqRecord.
    
    :param record: Biopython SeqRecord object.
    :param gene_identifiers: List of gene feature types to look for.
    :return: List of CDS SeqFeature objects.
    .. note:: antiSMASH gene features are usually called 'CDS'.
    """
    if gene_identifiers is None:
        gene_identifiers = ["CDS"]

    return [f for f in (record.features or []) if f.type in gene_identifiers]


def _iter_domains(record: "SeqRecord", domain_identifiers: list[str] | None = None) -> list["SeqFeature"]:
    """
    Iterate over domain features in a Biopython SeqRecord.
    
    :param record: Biopython SeqRecord object.
    :param domain_identifiers: List of domain feature types to look for.
    :return: List of domain SeqFeature objects.
    .. note:: antiSMASH domain features are usually called 'aSDomain'.
    """
    if domain_identifiers is None:
        domain_identifiers = ["aSDomain"]

    return [f for f in (record.features or []) if f.type in domain_identifiers]


def _start_end(feat: "SeqFeature") -> tuple[Strand, int, int]:
    """
    Return (strand, start, end) as a tuple from a SeqFeature.
    
    :param feat: Biopython SeqFeature object.
    :return: Tuple of (strand, start, end).
    :raises ValueError: If strand is not 1 or -1.
    .. note:: Location is 0-based, end-exclusive.
    """
    loc: "FeatureLocation" = feat.location

    match loc.strand:
        case 1:
            strand = Strand.FORWARD
        case -1:
            strand = Strand.REVERSE
        case _:
            raise ValueError(f"unexpected strand value: {loc.strand}")

    return strand, int(loc.start), int(loc.end)


def _in_bounds(child: "SeqFeature", parent: "SeqFeature") -> bool:
    """
    Check if a child feature is within the bounds of a parent feature.
    
    :param child: Child SeqFeature object.
    :param parent: Parent SeqFeature object.
    :return: True if child is within parent bounds, else False.
    """
    _, cs, ce = _start_end(child)
    _, ps, pe = _start_end(parent)

    return ps <= cs and ce <= pe


def _q1(feat: "SeqFeature", keys: Iterable[str]) -> str | None:
    """
    Return the first available qualifier value for any of the given keys, else None.

    :param feat: Biopython SeqFeature object.
    :param keys: Iterable of qualifier keys to check.
    :return: Qualifier value or None.
    """
    for k in keys:
        vals = feat.qualifiers.get(k)
        if vals:
            return vals[0]
    
    return None


def _gene_name(feat: "SeqFeature") -> str:
    """
    Get a gene name from a SeqFeature, or generate a default if none found.

    :param feat: Biopython SeqFeature object.
    :return: Gene name.
    """
    name = _q1(feat, ("locus_tag", "gene", "protein_id", "Name"))

    if name:
        return name
    
    strand, s, e = _start_end(feat)

    return f"CDS_{s}_{e}_{strand.value}"


def _gene_rec_from_feat(feat: "SeqFeature") -> Gene:
    """
    Create a Gene instance from a SeqFeature.

    :param feat: Biopython SeqFeature object.
    :return: Gene instance.
    :raises ValueError: If gene is not a Gene instance.
    """
    strand, s, e = _start_end(feat)
    name = _gene_name(feat)
    aa_seq = _q1(feat, ("translation",))

    if aa_seq is None:
        raise ValueError(f"Gene {name} has no sequence.")

    return Gene(
        id=name,
        start=s,
        end=e,
        strand=strand,
        sequence=aa_seq,
    )


def _domain_rec_from_feat(feat: "SeqFeature") -> Domain:
    """
    Create a Domain instance from a SeqFeature.

    :param feat: Biopython SeqFeature object.
    :return: Domain instance.
    :raises ValueError: If gene has no name, kind, or sequence.
    """
    _, s, e = _start_end(feat)
    name = _q1(feat, ("label", "product", "note"))
    kind = _q1(feat, ("aSDomain", "domain", "label"))
    aa_seq = _q1(feat, ("translation",))

    if name is None:
        raise ValueError(f"Gene {name} has no domain.")

    if kind is None:
        raise ValueError(f"Gene {name} has no domain.")

    if aa_seq is None:
        raise ValueError(f"Gene {name} has no domain.")

    return Domain(
        id=name,
        type=kind,
        start=s,
        end=e,
        sequence=aa_seq,
        raw_qualifiers={k: v for k, v in feat.qualifiers.items()},
    )


def collect_antismash_regions(
    input_file_path: Path | str,
    record: "SeqRecord",
    options: AntiSmashOptions
) -> list[Region]:
    """
    Collect antiSMASH regions from a GenBank record.

    :param input_file_path: Path to the GenBank file.
    :param record: GenBank record.
    :param options: Parsing options.
    :return: List of parsed regions.
    :raises ValueError: If record has no id.
    """
    # Quick check if record does contain id
    if record.id is None:
        raise ValueError(f"GenBank record has no ID.")

    # Get the file name without extension
    if isinstance(input_file_path, str):
        input_file_path = Path(input_file_path)
    file_name = Path(input_file_path).stem

    regions = _iter_regions(record, readout_level=options.readout_level)
    cds_list = _iter_cds(record, gene_identifiers=options.gene_identifiers)
    dom_list = _iter_domains(record, domain_identifiers=options.domain_identifiers)

    region_recs: list[Region] = []

    for reg in regions:
        _, rs, re = _start_end(reg)

        # Search for wanted qualifiers
        found_qualifiers: dict[str, Any] = {}
        for wanted in options.wanted_qualifiers or []:
            qualifiers = reg.qualifiers.get(wanted, []) or []
            found_qualifiers[wanted] = qualifiers

        # Genes inside region, sorted by coordinates
        gene_feats = [g for g in cds_list if _in_bounds(g, reg)]
        gene_feats.sort(key=lambda gf: (int(gf.location.start), int(gf.location.end)))

        # Make gene recs
        genes: list[Gene] = []
        for gf in gene_feats:
            g = _gene_rec_from_feat(gf)

            # Domains inside this gene, sorted by genomic start
            gene_doms = [df for df in dom_list if _in_bounds(df, gf)]
            gene_doms.sort(key=lambda df: (int(df.location.start), int(df.location.end)))
            dom_recs = [_domain_rec_from_feat(dd) for dd in gene_doms]

            g.domains = dom_recs
            genes.append(g)

        region_recs.append(Region(
            id=record.id,
            file_name=file_name,
            start=rs,
            end=re,
            qualifiers=found_qualifiers,
            genes=genes,
        ))

    return region_recs


def parse_antismash_gbk(path: Path, options: AntiSmashOptions) -> list:
    """
    Parse an antiSMASH-specific GenBank file.
    
    :param path: Path to the antiSMASH GenBank file.
    :param options: Parsing options.
    :return: List of parsed records.
    """
    log.info(f"parsing antiSMASH GenBank file: {path}")

    if not BIOPYTHON_AVAILABLE:
        raise ImportError("BioPython not available")

    out: list[Region] = []

    with open(path) as handle:
        for record in SeqIO.parse(handle, "genbank"):
            out.extend(collect_antismash_regions(path, record, options))

    return out
    

def load_regions(path: Path | str, options: RegionLoadOptions) -> list[Region]:
    """
    Load biosynthetic regions from a GenBank file.
    
    :param path: Path to the GenBank file.
    :param options: Options for loading regions, including source type and parsing parameters.
    :return: List of biosynthetic regions.
    :raises NotImplementedError: If the source is not implemented.
    """

    if isinstance(path, str):
        path = Path(path)

    match options:
        case AntiSmashOptions():
            return parse_antismash_gbk(Path(path), options)
        case _:
            raise NotImplementedError(f"loading regions from source '{options}' is not implemented")
