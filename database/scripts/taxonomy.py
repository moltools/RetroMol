"""NCBI taxonomy dump download + name/taxid resolution for phylogeny annotation.

Downloads and parses `taxdump.tar.gz` (https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz)
once, then resolves organism names <-> NCBI taxids and classifies a taxid's broad "type"
(Bacterium/Archaeon/Fungus/Other) by walking its lineage. MIBiG's JSON gives a taxid
directly (see extract_mibig_compounds.py's `taxonomy.ncbiTaxId`); NPAtlas gives only
genus/species text, resolved here by scientific-name/synonym lookup instead -- so both
sources end up with identically-standardized taxids and canonical names.
"""

from __future__ import annotations

import logging
import tarfile
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

from tqdm import tqdm

log = logging.getLogger(__name__)

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Fixed NCBI taxids for the lineage nodes used to classify a taxon's broad "type" (see
# TaxonomyDB.type_label_and_taxid) -- these are stable anchor points in NCBI's taxonomy,
# not resolved by name.
TAXID_BACTERIA = 2
TAXID_ARCHAEA = 2157
TAXID_FUNGI = 4751
TYPE_LABELS_BY_TAXID = {
    TAXID_BACTERIA: "Bacterium",
    TAXID_ARCHAEA: "Archaeon",
    TAXID_FUNGI: "Fungus",
}


def download_taxdump(dest_dir: str | Path, *, force: bool = False) -> Path:
    """Download and extract names.dmp/nodes.dmp from NCBI's taxdump into `dest_dir`
    (no-op if both files already exist there, unless `force`)."""
    dest_dir = Path(dest_dir).expanduser()
    dest_dir.mkdir(parents=True, exist_ok=True)

    names_path = dest_dir / "names.dmp"
    nodes_path = dest_dir / "nodes.dmp"
    if not force and names_path.exists() and nodes_path.exists():
        return dest_dir

    log.info("downloading NCBI taxdump from %s", TAXDUMP_URL)
    archive_path = dest_dir / "taxdump.tar.gz"
    with urllib.request.urlopen(TAXDUMP_URL) as resp, open(archive_path, "wb") as out:
        total = int(resp.headers.get("Content-Length") or 0) or None
        with tqdm(total=total, desc="download_taxdump", unit="B", unit_scale=True, unit_divisor=1024) as pbar:
            while chunk := resp.read(1024 * 1024):
                out.write(chunk)
                pbar.update(len(chunk))

    with tarfile.open(archive_path, mode="r:gz") as tar:
        for member in tar.getmembers():
            if member.name in ("names.dmp", "nodes.dmp"):
                tar.extract(member, path=dest_dir)
    archive_path.unlink()

    if not names_path.exists() or not nodes_path.exists():
        raise FileNotFoundError(f"taxdump at {TAXDUMP_URL} did not contain names.dmp/nodes.dmp")

    return dest_dir


def _iter_dmp_rows(path: Path) -> Iterator[list[str]]:
    """NCBI .dmp rows are "\\t|\\t"-separated, terminated by a trailing "\\t|"."""
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\t|")
            yield [field.strip() for field in line.split("\t|\t")]


@dataclass(frozen=True)
class PhylogenyResolution:
    type_label: str | None
    type_taxid: str | None
    genus: str | None
    genus_taxid: str | None
    species: str | None
    species_taxid: str | None


@dataclass
class TaxonomyDB:
    """In-memory index over NCBI's taxdump, built once per pipeline run."""

    parent_by_taxid: dict[int, int]
    rank_by_taxid: dict[int, str]
    scientific_name_by_taxid: dict[int, str]
    taxid_by_name: dict[str, int]  # lowercased scientific name/synonym -> taxid

    @classmethod
    def load(cls, taxdump_dir: str | Path) -> "TaxonomyDB":
        taxdump_dir = Path(taxdump_dir).expanduser()

        parent_by_taxid: dict[int, int] = {}
        rank_by_taxid: dict[int, str] = {}
        for row in _iter_dmp_rows(taxdump_dir / "nodes.dmp"):
            taxid, parent_taxid, rank = int(row[0]), int(row[1]), row[2]
            parent_by_taxid[taxid] = parent_taxid
            rank_by_taxid[taxid] = rank

        scientific_name_by_taxid: dict[int, str] = {}
        taxid_by_name: dict[str, int] = {}
        for row in _iter_dmp_rows(taxdump_dir / "names.dmp"):
            taxid, name_txt, name_class = int(row[0]), row[1], row[3]
            if name_class == "scientific name":
                scientific_name_by_taxid[taxid] = name_txt
            # Index every name class (scientific name, synonym, common name, ...) --
            # an organism string from a compound source can be any of these.
            taxid_by_name.setdefault(name_txt.lower(), taxid)

        log.info("loaded NCBI taxdump: %d nodes, %d names", len(parent_by_taxid), len(taxid_by_name))
        return cls(
            parent_by_taxid=parent_by_taxid,
            rank_by_taxid=rank_by_taxid,
            scientific_name_by_taxid=scientific_name_by_taxid,
            taxid_by_name=taxid_by_name,
        )

    def resolve_taxid(self, name: str | None) -> int | None:
        """Look up a taxid by scientific name or synonym, case-insensitive."""
        if not name:
            return None
        return self.taxid_by_name.get(name.strip().lower())

    def canonical_name(self, taxid: int | None) -> str | None:
        if taxid is None:
            return None
        return self.scientific_name_by_taxid.get(taxid)

    def lineage(self, taxid: int) -> list[int]:
        """The ancestor chain from `taxid` (inclusive) up to the root, stopping early
        if a cycle or a gap in nodes.dmp is hit."""
        chain = [taxid]
        seen = {taxid}
        current = taxid
        while current in self.parent_by_taxid and current != 1:
            parent = self.parent_by_taxid[current]
            if parent == current or parent in seen:
                break
            chain.append(parent)
            seen.add(parent)
            current = parent
        return chain

    def ancestor_at_rank(self, taxid: int | None, rank: str) -> int | None:
        """The ancestor of `taxid` (inclusive) whose rank is exactly `rank` (e.g.
        "genus", "species"), or None if no such ancestor exists in the lineage."""
        if taxid is None:
            return None
        for ancestor in self.lineage(taxid):
            if self.rank_by_taxid.get(ancestor) == rank:
                return ancestor
        return None

    def type_label_and_taxid(self, taxid: int | None) -> tuple[str | None, int | None]:
        """Classify a taxid's broad type by walking its lineage for one of the fixed
        Bacteria/Archaea/Fungi anchor nodes; "Other" if the lineage resolves but hits
        none of them (e.g. Viruses, Protista)."""
        if taxid is None:
            return None, None
        for ancestor in self.lineage(taxid):
            if ancestor in TYPE_LABELS_BY_TAXID:
                return TYPE_LABELS_BY_TAXID[ancestor], ancestor
        return "Other", None


def resolve_phylogeny(
    taxdb: TaxonomyDB | None,
    *,
    ncbi_tax_id: str | int | None = None,
    genus: str | None = None,
    species: str | None = None,
    fallback_type_label: str | None = None,
) -> PhylogenyResolution:
    """Standardize phylogeny fields to NCBI taxids where possible.

    If `ncbi_tax_id` is given (MIBiG), it's the source of truth: genus/species/type and
    their taxids are all derived from its lineage, overriding any given genus/species
    text. Otherwise (NPAtlas), `genus`/`species` text is resolved to a taxid by name
    lookup -- "Genus species" first, falling back to genus alone -- and, if resolved,
    canonical names/taxids for every rank are read back off the same lineage, so both
    sources end up identically standardized. If nothing resolves (no taxdb, or the
    name/taxid isn't found), the given genus/species/fallback_type_label pass through
    unchanged with no taxids -- unresolved is not an error here, just unenriched.
    """
    if taxdb is None:
        return PhylogenyResolution(fallback_type_label, None, genus, None, species, None)

    leaf_taxid: int | None = None
    if ncbi_tax_id:
        try:
            leaf_taxid = int(ncbi_tax_id)
        except (TypeError, ValueError):
            leaf_taxid = None
    elif genus:
        query = f"{genus} {species}" if species else genus
        leaf_taxid = taxdb.resolve_taxid(query) or taxdb.resolve_taxid(genus)

    if leaf_taxid is None:
        return PhylogenyResolution(fallback_type_label, None, genus, None, species, None)

    type_label, type_taxid = taxdb.type_label_and_taxid(leaf_taxid)

    genus_taxid = taxdb.ancestor_at_rank(leaf_taxid, "genus")
    genus_name = taxdb.canonical_name(genus_taxid) if genus_taxid else genus

    species_taxid = taxdb.ancestor_at_rank(leaf_taxid, "species")
    species_name = taxdb.canonical_name(species_taxid) if species_taxid else species

    return PhylogenyResolution(
        type_label=type_label or fallback_type_label,
        type_taxid=str(type_taxid) if type_taxid is not None else None,
        genus=genus_name,
        genus_taxid=str(genus_taxid) if genus_taxid is not None else None,
        species=species_name,
        species_taxid=str(species_taxid) if species_taxid is not None else None,
    )
