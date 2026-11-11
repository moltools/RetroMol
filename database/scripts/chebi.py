"""ChEBI role-ontology client for bioactivity annotation.

Downloads ChEBI's flat-file bulk release (ftp.ebi.ac.uk) once, then looks up a
compound's biological/chemical roles by its standard InChIKey. ChEBI's role ontology
(`has_role` edges under CHEBI:24432 "biological role" / CHEBI:51086 "chemical role")
is annotated on generic compound classes more often than on specific stereo-defined
structures -- confirmed
live (2026-08-26): erythromycin A (CHEBI:42355) itself carries no has_role edges, only
its `is_a` parent "erythromycin" (CHEBI:48923) does (roles: xenobiotic, bacterial
metabolite, environmental contaminant). So lookups walk a few steps up the `is_a`
ancestor chain collecting has_role edges at every level, not just the leaf's own.

Schema confirmed live by downloading and inspecting
ftp.ebi.ac.uk/pub/databases/chebi/flat_files/: structures.tsv (compound_id,
standard_inchi_key), relation.tsv (relation_type_id 4=has_role, 5=is_a; init_id/final_id
are compounds.tsv ids), compounds.tsv (id, name, chebi_accession).
"""

from __future__ import annotations

import csv
import gzip
import logging
import re
import urllib.request
from dataclasses import dataclass
from pathlib import Path

from tqdm import tqdm

log = logging.getLogger(__name__)

CHEBI_FLAT_FILES_URL = "https://ftp.ebi.ac.uk/pub/databases/chebi/flat_files"
CHEBI_FILES = ["compounds.tsv.gz", "structures.tsv.gz", "relation.tsv.gz"]

# ChEBI's `name` column carries inline HTML for italicized genus/species names and
# similar formatting (e.g. "<em>Aspergillus</em> metabolite", "<i>S</i>-configuration"
# -- confirmed live in compounds.tsv.gz/chemical_data.tsv.gz). Stripped here, once, at
# load time -- so every consumer of name_by_compound_id (role labels today, anything
# else later) gets plain text without needing to know this quirk exists.
_HTML_TAG_RE = re.compile(r"<[^>]+>")


def _strip_html(text: str) -> str:
    return _HTML_TAG_RE.sub("", text)

RELATION_TYPE_HAS_ROLE = "4"
RELATION_TYPE_IS_A = "5"

# ChEBI's role-ontology roots (CHEBI:50906 "role"'s direct children) -- a has_role
# target is classified by walking its own is_a ancestry for one of these two ids.
# "application" (CHEBI:33232), the third root, isn't surfaced -- not requested.
CHEBI_ID_BIOLOGICAL_ROLE = "24432"
CHEBI_ID_CHEMICAL_ROLE = "51086"

# How far to walk up `is_a` ancestry: role classification (shallow ontology, ~1-4 hops
# in practice) gets a bit more headroom than role-collection on the queried compound
# (kept shallow deliberately -- climbing further starts pulling in unrelated purely
# structural classes, e.g. erythromycin A's other is_a parent "cyclic ketone", which
# carry no roles of their own but would otherwise cost unbounded fan-out for nothing).
MAX_ROLE_CLASSIFICATION_HOPS = 10
MAX_ROLE_COLLECTION_HOPS = 4


def download_chebi_flat_files(dest_dir: str | Path, *, force: bool = False) -> Path:
    """Download compounds.tsv.gz/structures.tsv.gz/relation.tsv.gz into `dest_dir`
    (no-op if all three already exist, unless `force`)."""
    dest_dir = Path(dest_dir).expanduser()
    dest_dir.mkdir(parents=True, exist_ok=True)

    if not force and all((dest_dir / name).exists() for name in CHEBI_FILES):
        return dest_dir

    for name in CHEBI_FILES:
        url = f"{CHEBI_FLAT_FILES_URL}/{name}"
        dest = dest_dir / name
        log.info("downloading %s", url)
        with urllib.request.urlopen(url) as resp, open(dest, "wb") as out:
            total = int(resp.headers.get("Content-Length") or 0) or None
            with tqdm(total=total, desc=f"download_chebi[{name}]", unit="B", unit_scale=True, unit_divisor=1024) as pbar:
                while chunk := resp.read(1024 * 1024):
                    out.write(chunk)
                    pbar.update(len(chunk))

    return dest_dir


@dataclass(frozen=True)
class ChebiRole:
    label: str
    chebi_accession: str  # e.g. "CHEBI:35703" -- the role term's own id, for linking out


@dataclass(frozen=True)
class ChebiRoles:
    chebi_id: str  # the matched compound's own CHEBI accession, e.g. "CHEBI:42355"
    biological_roles: list[ChebiRole]
    chemical_roles: list[ChebiRole]


class ChebiDB:
    """In-memory index over ChEBI's flat files, built once per pipeline run."""

    def __init__(
        self,
        *,
        chebi_accession_by_compound_id: dict[str, str],
        name_by_compound_id: dict[str, str],
        compound_id_by_inchikey: dict[str, str],
        is_a_parents: dict[str, list[str]],
        has_role: dict[str, list[str]],
    ) -> None:
        self._chebi_accession_by_compound_id = chebi_accession_by_compound_id
        self._name_by_compound_id = name_by_compound_id
        self._compound_id_by_inchikey = compound_id_by_inchikey
        self._is_a_parents = is_a_parents
        self._has_role = has_role

    @classmethod
    def load(cls, chebi_dir: str | Path) -> "ChebiDB":
        chebi_dir = Path(chebi_dir).expanduser()

        chebi_accession_by_compound_id: dict[str, str] = {}
        name_by_compound_id: dict[str, str] = {}
        with gzip.open(chebi_dir / "compounds.tsv.gz", "rt", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh, delimiter="\t", quotechar='"')
            header = next(reader)
            id_idx = header.index("id")
            name_idx = header.index("name")
            accession_idx = header.index("chebi_accession")
            for row in reader:
                if len(row) <= accession_idx or not row[accession_idx]:
                    continue
                chebi_accession_by_compound_id[row[id_idx]] = row[accession_idx]
                name_by_compound_id[row[id_idx]] = _strip_html(row[name_idx])

        compound_id_by_inchikey: dict[str, str] = {}
        with gzip.open(chebi_dir / "structures.tsv.gz", "rt", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh, delimiter="\t", quotechar='"')
            header = next(reader)
            cid_idx = header.index("compound_id")
            key_idx = header.index("standard_inchi_key")
            for row in reader:
                if len(row) <= key_idx or not row[key_idx]:
                    continue
                compound_id_by_inchikey.setdefault(row[key_idx], row[cid_idx])

        is_a_parents: dict[str, list[str]] = {}
        has_role: dict[str, list[str]] = {}
        with gzip.open(chebi_dir / "relation.tsv.gz", "rt", encoding="utf-8", errors="replace") as fh:
            reader = csv.reader(fh, delimiter="\t", quotechar='"')
            header = next(reader)
            type_idx = header.index("relation_type_id")
            init_idx = header.index("init_id")
            final_idx = header.index("final_id")
            for row in reader:
                if len(row) <= final_idx:
                    continue
                rtype, init_id, final_id = row[type_idx], row[init_idx], row[final_idx]
                if rtype == RELATION_TYPE_IS_A:
                    is_a_parents.setdefault(init_id, []).append(final_id)
                elif rtype == RELATION_TYPE_HAS_ROLE:
                    has_role.setdefault(init_id, []).append(final_id)

        log.info(
            "loaded ChEBI flat files: %d compounds, %d structures, %d is_a edges, %d has_role edges",
            len(chebi_accession_by_compound_id), len(compound_id_by_inchikey), len(is_a_parents), len(has_role),
        )
        return cls(
            chebi_accession_by_compound_id=chebi_accession_by_compound_id,
            name_by_compound_id=name_by_compound_id,
            compound_id_by_inchikey=compound_id_by_inchikey,
            is_a_parents=is_a_parents,
            has_role=has_role,
        )

    def _walk_is_a(self, start: str, max_hops: int):
        """Yield ancestor compound ids reachable from `start` via `is_a`, breadth-first,
        up to `max_hops` levels, never revisiting a node."""
        seen = {start}
        frontier = [start]
        for _ in range(max_hops):
            next_frontier = []
            for cid in frontier:
                for parent in self._is_a_parents.get(cid, []):
                    if parent not in seen:
                        seen.add(parent)
                        next_frontier.append(parent)
                        yield parent
            if not next_frontier:
                return
            frontier = next_frontier

    def _role_category(self, role_compound_id: str) -> str | None:
        """Classify a has_role target by walking its own is_a ancestry for one of
        CHEBI_ID_BIOLOGICAL_ROLE/CHEBI_ID_CHEMICAL_ROLE. None if neither is hit (e.g.
        it falls under "application" instead, or the chain doesn't resolve)."""
        if role_compound_id == CHEBI_ID_BIOLOGICAL_ROLE:
            return "biological_role"
        if role_compound_id == CHEBI_ID_CHEMICAL_ROLE:
            return "chemical_role"
        for ancestor in self._walk_is_a(role_compound_id, MAX_ROLE_CLASSIFICATION_HOPS):
            if ancestor == CHEBI_ID_BIOLOGICAL_ROLE:
                return "biological_role"
            if ancestor == CHEBI_ID_CHEMICAL_ROLE:
                return "chemical_role"
        return None

    def roles_for_inchikey(self, inchikey: str) -> ChebiRoles | None:
        """None if `inchikey` has no match in ChEBI at all. Roles are collected from
        the matched compound and a few levels of its `is_a` ancestors (see module
        docstring -- ChEBI usually annotates roles on a generic parent class, not every
        specific stereo-defined child structure)."""
        compound_id = self._compound_id_by_inchikey.get(inchikey)
        if compound_id is None:
            return None

        role_compound_ids: set[str] = set(self._has_role.get(compound_id, []))
        for ancestor in self._walk_is_a(compound_id, MAX_ROLE_COLLECTION_HOPS):
            role_compound_ids.update(self._has_role.get(ancestor, []))

        biological_roles: list[ChebiRole] = []
        chemical_roles: list[ChebiRole] = []
        for role_id in role_compound_ids:
            label = self._name_by_compound_id.get(role_id)
            accession = self._chebi_accession_by_compound_id.get(role_id)
            if not label or not accession:
                continue
            category = self._role_category(role_id)
            if category == "biological_role":
                biological_roles.append(ChebiRole(label=label, chebi_accession=accession))
            elif category == "chemical_role":
                chemical_roles.append(ChebiRole(label=label, chebi_accession=accession))

        return ChebiRoles(
            chebi_id=self._chebi_accession_by_compound_id.get(compound_id, compound_id),
            biological_roles=biological_roles,
            chemical_roles=chemical_roles,
        )
