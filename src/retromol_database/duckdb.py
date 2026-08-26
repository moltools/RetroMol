from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Literal, Sequence

import duckdb
import numpy as np

from retromol_fingerprint.fingerprint import TOKEN_LINK, TOKEN_UNK

ENTRY_TYPES = ("compound", "bgc")
EntryType = Literal["compound", "bgc"]
FINGERPRINT_SIZE = 1024


@dataclass(frozen=True)
class EntrySource:
    name: str
    database_name: str
    url: str | None


@dataclass(frozen=True)
class Entry:
    id: str
    # Primary display name/url -- the first-ever-inserted source for this entry (see
    # `sources` below and entry_sources.seq). Kept as plain fields (not derived via a
    # property) so every existing read site and the synthetic upload-candidate Entry(...)
    # constructions in routes/discovery.py keep working unchanged.
    name: str
    url: str | None
    raw: str | None
    type: EntryType
    primary_sequence: list[str]
    fingerprint: list[float]
    # Every (name, database_name, url) a source has contributed for this entry, ordered
    # by insertion. Empty for synthetic (non-database) entries, e.g. session uploads.
    sources: list[EntrySource] = field(default_factory=list)


@dataclass(frozen=True)
class SearchResult:
    entry: Entry
    similarity: float


@dataclass(frozen=True)
class Count:
    label: str
    count: int


@dataclass(frozen=True)
class DatabaseStats:
    total_entries: int
    counts_by_type: list[Count]
    sequence_length_min: int
    sequence_length_max: int
    sequence_length_avg: float
    unique_block_count: int
    with_source_url_count: int
    without_source_url_count: int


@dataclass(frozen=True)
class AnnotationTerm:
    id: str
    category: str
    rank: str | None
    label: str
    parent_id: str | None
    # An id in whatever external database this term comes from (NCBI taxid for
    # phylogeny; ATC code/ChEBI accession/ChEMBL id for bioactivity), for linking out
    # to that database's own page -- None for categories with no such external page
    # (biosynthetic_class, chemical_class).
    external_id: str | None = None


@dataclass(frozen=True)
class AnnotationStats:
    with_annotation_count: int
    without_annotation_count: int
    counts_by_category: list[Count]
    # Phylogeny: one chart per rank.
    phylogeny_type_counts: list[Count]
    phylogeny_genus_counts: list[Count]
    phylogeny_species_counts: list[Count]
    # MIBiG's own coarse biosynthetic-class label -- a distinct category from chemical_class
    # below, not a rank within it (see biosynthetic_class_annotations table comment).
    biosynthetic_class_counts: list[Count]
    # Chemical class: NPClassifier's three structure-derived levels (see
    # database/scripts/annotate_npclassifier.py), kept as separate charts rather than
    # pooled, since they're different ranks within the same classification scheme.
    chemical_class_pathway_counts: list[Count]
    chemical_class_superclass_counts: list[Count]
    chemical_class_class_counts: list[Count]
    # Bioactivity: ChEMBL's two signals (see database/scripts/annotate_chembl.py) --
    # WHO ATC therapeutic category and clinical-development phase -- plus ChEBI's role
    # ontology (see database/scripts/annotate_chebi.py).
    bioactivity_atc_counts: list[Count]
    bioactivity_max_phase_counts: list[Count]
    bioactivity_biological_role_counts: list[Count]
    bioactivity_chemical_role_counts: list[Count]


def _normalize_entry_type(entry_type: str) -> EntryType:
    if entry_type not in ENTRY_TYPES:
        raise ValueError(f"entry type must be one of {ENTRY_TYPES}, got {entry_type}")

    return entry_type  # type: ignore[return-value]


def _normalize_primary_sequence(primary_sequence: Sequence[str]) -> list[str]:
    return [str(token) for token in primary_sequence]


def _normalize_fingerprint(fingerprint: Sequence[float] | np.ndarray) -> list[float]:
    fp = np.asarray(fingerprint, dtype=np.float32)

    if fp.shape != (FINGERPRINT_SIZE,):
        raise ValueError(f"fingerprint must have shape ({FINGERPRINT_SIZE},), got {fp.shape}")

    return fp.astype(float).tolist()


class RetroMolDuckDB:
    def __init__(self, path: str | Path, *, read_only: bool = False) -> None:
        self.path = Path(path).expanduser()
        self.con = duckdb.connect(str(self.path), read_only=read_only)

    def __enter__(self) -> "RetroMolDuckDB":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()

    @classmethod
    def create(cls, path: str | Path, *, overwrite: bool = False) -> "RetroMolDuckDB":
        path: Path = Path(path).expanduser()

        if overwrite and path.exists():
            path.unlink()

        db = cls(path)
        db.create_schema()
        return db

    @classmethod
    def open(cls, path: str | Path, *, read_only: bool = False) -> "RetroMolDuckDB":
        path: Path = Path(path).expanduser()

        if read_only and not path.exists():
            raise FileNotFoundError(path)

        db = cls(path, read_only=read_only)

        if not read_only:
            db.create_schema()

        return db

    def close(self) -> None:
        self.con.close()

    def create_schema(self) -> None:
        self.con.execute("CREATE SEQUENCE IF NOT EXISTS entry_sources_seq")
        self.con.execute(
            f"""
            CREATE TABLE IF NOT EXISTS entries (
                id VARCHAR PRIMARY KEY,
                type VARCHAR NOT NULL CHECK (type IN ('compound', 'bgc')),
                raw VARCHAR,
                content_hash VARCHAR,
                primary_sequence VARCHAR[] NOT NULL,
                fingerprint FLOAT[{FINGERPRINT_SIZE}] NOT NULL
            )
            """
        )
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS entry_sources (
                seq BIGINT DEFAULT nextval('entry_sources_seq'),
                entry_id VARCHAR NOT NULL,
                name VARCHAR NOT NULL,
                database_name VARCHAR NOT NULL,
                url VARCHAR,
                PRIMARY KEY (entry_id, database_name, name)
            )
            """
        )
        # Four dedicated annotation tables, one per category, each with its own typed
        # columns instead of a single generic label -- see RetroMolDuckDB.add_phylogeny_annotation
        # / add_bioactivity_annotation / add_biosynthetic_class_annotation / add_chemical_class_annotation.
        # Phylogeny is one row per entry (an entry has exactly one organism); the other three
        # are multi-valued per entry. biosynthetic_class/chemical_class are populated for
        # compounds only (they describe the molecule, not the organism/cluster) -- enforced
        # by pipeline callers, not by this schema. Bioactivity population logic isn't wired
        # up yet; that table exists ahead of that.
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS phylogeny_annotations (
                entry_id VARCHAR PRIMARY KEY,
                type_label VARCHAR,
                type_taxid VARCHAR,
                genus VARCHAR,
                genus_taxid VARCHAR,
                species VARCHAR,
                species_taxid VARCHAR
            )
            """
        )
        # `level` distinguishes bioactivity's four signals (see
        # database/scripts/annotate_chembl.py / annotate_chebi.py): ChEMBL's
        # "chembl_max_phase" (clinical-development stage -- "Approved"/"Phase 3"/... --
        # only stored when meaningful, same presence-only convention as chemical_class's
        # is_glycoside) and "chembl_atc" (WHO ATC therapeutic category), plus ChEBI's
        # "chebi_biological_role" and "chebi_chemical_role" (its role ontology, under
        # CHEBI:24432/CHEBI:51086). All four are multi-valued per entry except
        # chembl_max_phase. `external_id` is the id in that row's own external database
        # (an ATC code, a ChEBI accession, or a ChEMBL id for the max_phase row's own
        # compound page) -- for building a "view on <source>" link; null when unresolved.
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS bioactivity_annotations (
                entry_id VARCHAR NOT NULL,
                level VARCHAR NOT NULL,
                label VARCHAR NOT NULL,
                external_id VARCHAR,
                PRIMARY KEY (entry_id, level, label)
            )
            """
        )
        # MIBiG's own coarse biosynthetic-class label (PKS/NRPS/RiPP/terpene/saccharide/other,
        # from the BGC's annotated biosynthesis machinery) -- a distinct classification scheme
        # from chemical_class below (NPClassifier's structure-derived classification of the
        # compound itself), not a rank/level within it.
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS biosynthetic_class_annotations (
                entry_id VARCHAR NOT NULL,
                label VARCHAR NOT NULL,
                PRIMARY KEY (entry_id, label)
            )
            """
        )
        # `level` distinguishes NPClassifier's "pathway"/"superclass"/"class"/"is_glycoside"
        # (structure-derived, from a compound's own SMILES -- see
        # database/scripts/annotate_npclassifier.py) -- every chemical_class row comes from
        # NPClassifier now (MIBiG's biosynthetic class lives in biosynthetic_class_annotations
        # above instead). The three list-valued levels can each carry more than one label per
        # entry (hence label being part of the primary key, same as bioactivity); is_glycoside
        # is stored as a single ("is_glycoside", "Yes") row only when true -- absence means
        # false/unclassified, the same presence-only convention every other annotation table uses.
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS chemical_class_annotations (
                entry_id VARCHAR NOT NULL,
                level VARCHAR NOT NULL,
                label VARCHAR NOT NULL,
                PRIMARY KEY (entry_id, level, label)
            )
            """
        )
        # Read-only views reconstructing the old generic (id/category/rank/label/parent_id)
        # term shape and (entry_id, term_id) links across all three tables above -- so
        # annotation_stats/enrichment queries (which are category-agnostic) don't need to
        # know about the dedicated tables. term_id encodings match the pre-migration scheme
        # (e.g. "phylogeny:genus:bacterium:streptomyces") so existing consumers are unaffected.
        self.con.execute(
            """
            CREATE OR REPLACE VIEW annotation_terms AS
            SELECT id, category, rank, label, parent_id, any_value(taxid) AS external_id
            FROM (
                SELECT
                    'phylogeny:type:' || lower(type_label) AS id,
                    'phylogeny' AS category, 'type' AS rank, type_label AS label,
                    CAST(NULL AS VARCHAR) AS parent_id, type_taxid AS taxid
                FROM phylogeny_annotations WHERE type_label IS NOT NULL
                UNION ALL
                SELECT
                    'phylogeny:genus:' || lower(type_label) || ':' || lower(genus) AS id,
                    'phylogeny' AS category, 'genus' AS rank, genus AS label,
                    'phylogeny:type:' || lower(type_label) AS parent_id, genus_taxid AS taxid
                FROM phylogeny_annotations WHERE type_label IS NOT NULL AND genus IS NOT NULL
                UNION ALL
                SELECT
                    'phylogeny:species:' || lower(type_label) || ':' || lower(genus) || ':' || lower(species) AS id,
                    'phylogeny' AS category, 'species' AS rank, species AS label,
                    'phylogeny:genus:' || lower(type_label) || ':' || lower(genus) AS parent_id, species_taxid AS taxid
                FROM phylogeny_annotations
                WHERE type_label IS NOT NULL AND genus IS NOT NULL AND species IS NOT NULL
                UNION ALL
                SELECT
                    'bioactivity:' || level || ':' || lower(label) AS id,
                    'bioactivity' AS category, level AS rank, label,
                    CAST(NULL AS VARCHAR) AS parent_id, external_id AS taxid
                FROM bioactivity_annotations
                UNION ALL
                SELECT
                    'biosynthetic_class:' || lower(label) AS id,
                    'biosynthetic_class' AS category, CAST(NULL AS VARCHAR) AS rank, label,
                    CAST(NULL AS VARCHAR) AS parent_id, CAST(NULL AS VARCHAR) AS taxid
                FROM biosynthetic_class_annotations
                UNION ALL
                SELECT
                    'chemical_class:' || level || ':' || lower(label) AS id,
                    'chemical_class' AS category, level AS rank, label,
                    CAST(NULL AS VARCHAR) AS parent_id, CAST(NULL AS VARCHAR) AS taxid
                FROM chemical_class_annotations
            )
            GROUP BY id, category, rank, label, parent_id
            """
        )
        self.con.execute(
            """
            CREATE OR REPLACE VIEW entry_annotations AS
            SELECT entry_id, 'phylogeny:type:' || lower(type_label) AS term_id
            FROM phylogeny_annotations WHERE type_label IS NOT NULL
            UNION ALL
            SELECT entry_id, 'phylogeny:genus:' || lower(type_label) || ':' || lower(genus)
            FROM phylogeny_annotations WHERE type_label IS NOT NULL AND genus IS NOT NULL
            UNION ALL
            SELECT entry_id, 'phylogeny:species:' || lower(type_label) || ':' || lower(genus) || ':' || lower(species)
            FROM phylogeny_annotations
            WHERE type_label IS NOT NULL AND genus IS NOT NULL AND species IS NOT NULL
            UNION ALL
            SELECT entry_id, 'bioactivity:' || level || ':' || lower(label) FROM bioactivity_annotations
            UNION ALL
            SELECT entry_id, 'biosynthetic_class:' || lower(label) FROM biosynthetic_class_annotations
            UNION ALL
            SELECT entry_id, 'chemical_class:' || level || ':' || lower(label) FROM chemical_class_annotations
            """
        )

    def add_entry(
        self,
        *,
        entry_id: str,
        name: str,
        database_name: str,
        url: str | None,
        raw: str | None,
        entry_type: str,
        primary_sequence: Sequence[str],
        fingerprint: Sequence[float] | np.ndarray,
        content_hash: str | None = None,
    ) -> str:
        """
        Add (or extend) an entry.

        `entry_id` is the caller-supplied molecular identity -- an InChIKey for a
        compound, or a hash of the source .gbk file's content plus region id for a
        bgc (see database/scripts/load_compounds.py / load_bgcs.py). If an entry
        with this id already exists, its stored `raw`/`primary_sequence`/
        `fingerprint` are left untouched (they describe the same molecule/BGC
        either way) and only a new `(name, database_name, url)` source row is
        added -- or, if that exact (entry_id, database_name, name) combination was
        already recorded, its url is refreshed.
        """
        entry_type = _normalize_entry_type(entry_type)
        sequence = _normalize_primary_sequence(primary_sequence)
        fp = _normalize_fingerprint(fingerprint)

        self.con.execute(
            """
            INSERT INTO entries (id, type, raw, content_hash, primary_sequence, fingerprint)
            VALUES (?, ?, ?, ?, ?, ?)
            ON CONFLICT (id) DO NOTHING
            """,
            [entry_id, entry_type, raw, content_hash, sequence, fp],
        )

        self.con.execute(
            """
            INSERT INTO entry_sources (entry_id, name, database_name, url)
            VALUES (?, ?, ?, ?)
            ON CONFLICT (entry_id, database_name, name) DO UPDATE SET url = excluded.url
            """,
            [entry_id, name, database_name, url],
        )

        return entry_id

    def bgc_content_hash_exists(self, content_hash: str) -> bool:
        """Whether a bgc entry sourced from a .gbk file with this content hash already exists."""
        row = self.con.execute(
            "SELECT 1 FROM entries WHERE type = 'bgc' AND content_hash = ? LIMIT 1",
            [content_hash],
        ).fetchone()
        return row is not None

    def count(self) -> int:
        return int(self.con.execute("SELECT count(*) FROM entries").fetchone()[0])

    def add_phylogeny_annotation(
        self,
        entry_id: str,
        *,
        type_label: str | None,
        type_taxid: str | None = None,
        genus: str | None = None,
        genus_taxid: str | None = None,
        species: str | None = None,
        species_taxid: str | None = None,
    ) -> None:
        """Set (or replace) `entry_id`'s single phylogeny row. `species`/`species_taxid`
        are only meaningful -- and only surfaced by the annotation_terms/entry_annotations
        views -- when `genus` is also given; a taxid is optional at every rank (e.g. NPAtlas
        names that don't resolve against the NCBI taxdump still store their raw genus/species
        text with a NULL taxid). No-ops if `type_label` isn't resolvable."""
        if not type_label:
            return

        self.con.execute(
            """
            INSERT INTO phylogeny_annotations
                (entry_id, type_label, type_taxid, genus, genus_taxid, species, species_taxid)
            VALUES (?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT (entry_id) DO UPDATE SET
                type_label = excluded.type_label,
                type_taxid = excluded.type_taxid,
                genus = excluded.genus,
                genus_taxid = excluded.genus_taxid,
                species = excluded.species,
                species_taxid = excluded.species_taxid
            """,
            [entry_id, type_label, type_taxid, genus, genus_taxid, species, species_taxid],
        )

    def add_bioactivity_annotation(
        self, entry_id: str, *, level: str, label: str, external_id: str | None = None
    ) -> None:
        """Link `entry_id` (a compound) to a bioactivity label at the given `level` (e.g.
        "chembl_max_phase", "chembl_atc", "chebi_biological_role", "chebi_chemical_role" --
        see database/scripts/annotate_chembl.py / annotate_chebi.py). Compounds only --
        callers must not use this for bgc entries."""
        self.con.execute(
            """
            INSERT INTO bioactivity_annotations (entry_id, level, label, external_id)
            VALUES (?, ?, ?, ?)
            ON CONFLICT (entry_id, level, label) DO NOTHING
            """,
            [entry_id, level, label, external_id],
        )

    def add_biosynthetic_class_annotation(self, entry_id: str, label: str) -> None:
        """Link `entry_id` (a compound) to a MIBiG biosynthetic-class label (PKS/NRPS/RiPP/...).
        Compounds only -- callers must not use this for bgc entries."""
        self.con.execute(
            """
            INSERT INTO biosynthetic_class_annotations (entry_id, label)
            VALUES (?, ?)
            ON CONFLICT (entry_id, label) DO NOTHING
            """,
            [entry_id, label],
        )

    def add_chemical_class_annotation(self, entry_id: str, *, level: str, label: str) -> None:
        """Link `entry_id` (a compound) to a chemical-class label at the given `level`
        (e.g. "biosyn_class", or NPClassifier's "pathway"/"superclass"/"class"/"is_glycoside" --
        see the chemical_class_annotations table comment in create_schema). Compounds only --
        callers must not use this for bgc entries."""
        self.con.execute(
            """
            INSERT INTO chemical_class_annotations (entry_id, level, label)
            VALUES (?, ?, ?)
            ON CONFLICT (entry_id, level, label) DO NOTHING
            """,
            [entry_id, level, label],
        )

    def count_entries_by_type(self, entry_types: Sequence[str]) -> int:
        if not entry_types:
            return 0
        types = [_normalize_entry_type(t) for t in entry_types]
        return int(
            self.con.execute(
                "SELECT count(*) FROM entries WHERE type IN (SELECT UNNEST(?))",
                [types],
            ).fetchone()[0]
        )

    def annotation_term_counts(self, entry_ids: Sequence[str]) -> dict[str, int]:
        """For every term linked to at least one of the given entry ids, how many of those
        ids carry it. Used as the "selected" side of an enrichment contingency table."""
        if not entry_ids:
            return {}

        rows = self.con.execute(
            """
            SELECT term_id, count(DISTINCT entry_id)
            FROM entry_annotations
            WHERE entry_id IN (SELECT UNNEST(?))
            GROUP BY term_id
            """,
            [list(entry_ids)],
        ).fetchall()
        return {str(row[0]): int(row[1]) for row in rows}

    def annotation_term_counts_for_types(self, entry_types: Sequence[str]) -> dict[str, int]:
        """For every term, how many entries of the given type(s) carry it -- the
        "background pool" side of an enrichment contingency table."""
        if not entry_types:
            return {}

        types = [_normalize_entry_type(t) for t in entry_types]
        rows = self.con.execute(
            """
            SELECT ea.term_id, count(DISTINCT ea.entry_id)
            FROM entry_annotations ea
            JOIN entries e ON e.id = ea.entry_id
            WHERE e.type IN (SELECT UNNEST(?))
            GROUP BY ea.term_id
            """,
            [types],
        ).fetchall()
        return {str(row[0]): int(row[1]) for row in rows}

    def annotation_terms_by_ids(self, term_ids: Sequence[str]) -> dict[str, AnnotationTerm]:
        if not term_ids:
            return {}

        rows = self.con.execute(
            """
            SELECT id, category, rank, label, parent_id, external_id
            FROM annotation_terms
            WHERE id IN (SELECT UNNEST(?))
            """,
            [list(term_ids)],
        ).fetchall()
        return {
            str(row[0]): AnnotationTerm(
                id=str(row[0]), category=str(row[1]), rank=row[2], label=str(row[3]), parent_id=row[4],
                external_id=row[5],
            )
            for row in rows
        }

    def entry_annotation_terms(self, entry_id: str) -> list[AnnotationTerm]:
        """Every annotation term linked to `entry_id`, across all four categories --
        the "show me everything known about this compound/bgc" query used by the
        Discovery tab's expanded result view."""
        rows = self.con.execute(
            """
            SELECT t.id, t.category, t.rank, t.label, t.parent_id, t.external_id
            FROM entry_annotations ea
            JOIN annotation_terms t ON t.id = ea.term_id
            WHERE ea.entry_id = ?
            ORDER BY t.category, t.rank NULLS FIRST, t.label
            """,
            [entry_id],
        ).fetchall()
        return [
            AnnotationTerm(
                id=str(row[0]), category=str(row[1]), rank=row[2], label=str(row[3]), parent_id=row[4],
                external_id=row[5],
            )
            for row in rows
        ]

    def search_entries(
        self, query: str, *, entry_type: str | None = None, limit: int = 100
    ) -> list[Entry]:
        """Look up up to `limit` entries by a case-insensitive substring match on any of
        their source names, or an exact match on their id -- the "query and select" lookup
        used by the Enrichment tab (distinct from `closest()`'s fingerprint-similarity search)."""
        if limit < 1:
            raise ValueError("limit must be >= 1")

        where_sql = "WHERE (e.id = $2 OR es.name ILIKE '%' || $2 || '%')"
        params: list[object] = [limit, query]

        if entry_type is not None:
            entry_type = _normalize_entry_type(entry_type)
            where_sql += " AND e.type = $3"
            params.append(entry_type)

        rows = self.con.execute(
            f"""
            SELECT DISTINCT e.id, e.raw, e.type, e.primary_sequence, e.fingerprint
            FROM entries e
            LEFT JOIN entry_sources es ON es.entry_id = e.id
            {where_sql}
            ORDER BY e.id
            LIMIT $1
            """,
            params,
        ).fetchall()

        sources_by_id = self._sources_for_entry_ids([str(row[0]) for row in rows])
        return [_entry_from_row(row, sources_by_id.get(str(row[0]), [])) for row in rows]

    def _label_counts(self, *, category: str, rank: str | None = None, limit: int | None = None) -> list[Count]:
        """Distinct-entry counts per label, for one (category, rank) slice of
        annotation_terms/entry_annotations -- the building block for every chart on the
        Dashboard's Annotations section."""
        where_sql = "WHERE t.category = ?"
        params: list[object] = [category]
        if rank is not None:
            where_sql += " AND t.rank = ?"
            params.append(rank)

        limit_sql = " LIMIT ?" if limit is not None else ""
        if limit is not None:
            params.append(limit)

        rows = self.con.execute(
            f"""
            SELECT t.label, count(DISTINCT ea.entry_id)
            FROM entry_annotations ea
            JOIN annotation_terms t ON t.id = ea.term_id
            {where_sql}
            GROUP BY t.label
            ORDER BY count(DISTINCT ea.entry_id) DESC
            {limit_sql}
            """,
            params,
        ).fetchall()
        return [Count(label=str(row[0]), count=int(row[1])) for row in rows]

    def annotation_stats(self) -> AnnotationStats:
        """Summary counts over annotation_terms/entry_annotations, for the Dashboard."""
        with_annotation_count = int(
            self.con.execute("SELECT count(DISTINCT entry_id) FROM entry_annotations").fetchone()[0]
        )
        without_annotation_count = self.count() - with_annotation_count

        category_rows = self.con.execute(
            """
            SELECT t.category, count(DISTINCT ea.entry_id)
            FROM entry_annotations ea
            JOIN annotation_terms t ON t.id = ea.term_id
            GROUP BY t.category
            ORDER BY t.category
            """
        ).fetchall()
        counts_by_category = [Count(label=str(row[0]), count=int(row[1])) for row in category_rows]

        return AnnotationStats(
            with_annotation_count=with_annotation_count,
            without_annotation_count=without_annotation_count,
            counts_by_category=counts_by_category,
            phylogeny_type_counts=self._label_counts(category="phylogeny", rank="type"),
            phylogeny_genus_counts=self._label_counts(category="phylogeny", rank="genus", limit=15),
            phylogeny_species_counts=self._label_counts(category="phylogeny", rank="species", limit=15),
            biosynthetic_class_counts=self._label_counts(category="biosynthetic_class"),
            chemical_class_pathway_counts=self._label_counts(category="chemical_class", rank="pathway", limit=15),
            chemical_class_superclass_counts=self._label_counts(category="chemical_class", rank="superclass", limit=15),
            chemical_class_class_counts=self._label_counts(category="chemical_class", rank="class", limit=15),
            bioactivity_atc_counts=self._label_counts(category="bioactivity", rank="chembl_atc", limit=15),
            bioactivity_max_phase_counts=self._label_counts(category="bioactivity", rank="chembl_max_phase"),
            bioactivity_biological_role_counts=self._label_counts(
                category="bioactivity", rank="chebi_biological_role", limit=15
            ),
            bioactivity_chemical_role_counts=self._label_counts(
                category="bioactivity", rank="chebi_chemical_role", limit=15
            ),
        )

    def stats(self) -> DatabaseStats:
        """
        Compute summary statistics over the whole entries table, for display on a
        dashboard/overview page.

        Building blocks are the tokens in each entry's primary_sequence -- e.g. amino
        acid names, PK reduction-state groups, or tailoring events like "methylation".
        TOKEN_UNK ("<UNK>") marks a block RetroMol couldn't identify and TOKEN_LINK
        ("<LINK>") just joins two merged paths within one entry's sequence -- neither
        is a real building block, so both are excluded from unique_block_count.

        :return: a DatabaseStats snapshot
        """
        total_entries = self.count()

        type_rows = self.con.execute(
            "SELECT type, count(*) FROM entries GROUP BY type ORDER BY type"
        ).fetchall()
        counts_by_type = [Count(label=str(row[0]), count=int(row[1])) for row in type_rows]

        length_row = self.con.execute(
            """
            SELECT
                min(len(primary_sequence)),
                max(len(primary_sequence)),
                avg(len(primary_sequence))
            FROM entries
            """
        ).fetchone()
        sequence_length_min = int(length_row[0]) if length_row[0] is not None else 0
        sequence_length_max = int(length_row[1]) if length_row[1] is not None else 0
        sequence_length_avg = float(length_row[2]) if length_row[2] is not None else 0.0

        unique_block_count = int(
            self.con.execute(
                """
                SELECT count(DISTINCT token)
                FROM (SELECT unnest(primary_sequence) AS token FROM entries)
                WHERE token NOT IN (?, ?)
                """,
                [TOKEN_UNK, TOKEN_LINK],
            ).fetchone()[0]
        )

        url_row = self.con.execute(
            """
            SELECT
                count(*) FILTER (WHERE has_url),
                count(*) FILTER (WHERE NOT has_url)
            FROM (
                SELECT e.id, bool_or(es.url IS NOT NULL) AS has_url
                FROM entries e
                LEFT JOIN entry_sources es ON es.entry_id = e.id
                GROUP BY e.id
            )
            """
        ).fetchone()
        with_source_url_count = int(url_row[0])
        without_source_url_count = int(url_row[1])

        return DatabaseStats(
            total_entries=total_entries,
            counts_by_type=counts_by_type,
            sequence_length_min=sequence_length_min,
            sequence_length_max=sequence_length_max,
            sequence_length_avg=sequence_length_avg,
            unique_block_count=unique_block_count,
            with_source_url_count=with_source_url_count,
            without_source_url_count=without_source_url_count,
        )

    def _sources_for_entry_ids(self, entry_ids: Sequence[str]) -> dict[str, list[EntrySource]]:
        """Batch-fetch every (name, database_name, url) source for the given entry ids, ordered by insertion."""
        if not entry_ids:
            return {}

        rows = self.con.execute(
            """
            SELECT entry_id, name, database_name, url
            FROM entry_sources
            WHERE entry_id IN (SELECT UNNEST(?))
            ORDER BY entry_id, seq
            """,
            [list(entry_ids)],
        ).fetchall()

        out: dict[str, list[EntrySource]] = {}
        for entry_id, name, database_name, url in rows:
            out.setdefault(str(entry_id), []).append(
                EntrySource(name=str(name), database_name=str(database_name), url=url)
            )
        return out

    def get_entry(self, entry_id: str) -> Entry | None:
        row = self.con.execute(
            """
            SELECT id, raw, type, primary_sequence, fingerprint
            FROM entries
            WHERE id = ?
            """,
            [entry_id],
        ).fetchone()

        if row is None:
            return None

        sources = self._sources_for_entry_ids([entry_id]).get(entry_id, [])
        return _entry_from_row(row, sources)

    def get_entries(self, entry_ids: Sequence[str]) -> list[Entry]:
        """Batch version of get_entry -- one round trip for up to `len(entry_ids)` entries,
        in no particular order (missing ids are silently omitted, not errored)."""
        if not entry_ids:
            return []

        rows = self.con.execute(
            """
            SELECT id, raw, type, primary_sequence, fingerprint
            FROM entries
            WHERE id IN (SELECT UNNEST(?))
            """,
            [list(entry_ids)],
        ).fetchall()

        sources_by_id = self._sources_for_entry_ids([str(row[0]) for row in rows])
        return [_entry_from_row(row, sources_by_id.get(str(row[0]), [])) for row in rows]

    def iter_entries(self) -> Iterator[Entry]:
        rows = self.con.execute(
            """
            SELECT id, raw, type, primary_sequence, fingerprint
            FROM entries
            ORDER BY id
            """
        ).fetchall()

        sources_by_id = self._sources_for_entry_ids([str(row[0]) for row in rows])
        for row in rows:
            yield _entry_from_row(row, sources_by_id.get(str(row[0]), []))

    def closest(
        self,
        fingerprint: Sequence[float] | np.ndarray,
        *,
        limit: int = 1000,
        entry_type: str | None = None,
    ) -> list[SearchResult]:
        if limit < 1:
            raise ValueError("limit must be >= 1")

        fp = _normalize_fingerprint(fingerprint)
        params: list[object] = [fp]
        where_sql = ""

        if entry_type is not None:
            entry_type: str = _normalize_entry_type(entry_type)
            where_sql = "WHERE type = ?"
            params.append(entry_type)

        params.append(limit)

        rows = self.con.execute(
            f"""
            SELECT
                id,
                raw,
                type,
                primary_sequence,
                fingerprint,
                array_cosine_similarity(fingerprint, ?::FLOAT[{FINGERPRINT_SIZE}]) AS similarity
            FROM entries
            {where_sql}
            ORDER BY similarity DESC
            LIMIT ?
            """,
            params,
        ).fetchall()

        sources_by_id = self._sources_for_entry_ids([str(row[0]) for row in rows])

        return [
            SearchResult(
                entry=_entry_from_row(row[:5], sources_by_id.get(str(row[0]), [])),
                similarity=float(row[5]),
            )
            for row in rows
        ]

    def export_parquet(self, path: str | Path) -> None:
        self.con.execute(
            "COPY entries TO ? (FORMAT PARQUET)",
            [str(path)],
        )

def _entry_from_row(row, sources: list[EntrySource]) -> Entry:
    primary = sources[0] if sources else None
    return Entry(
        id=str(row[0]),
        name=primary.name if primary else str(row[0]),
        url=primary.url if primary else None,
        raw=row[1],
        type=_normalize_entry_type(row[2]),
        primary_sequence=list(row[3]),
        fingerprint=list(row[4]),
        sources=sources,
    )


def create_database(path: str | Path, *, overwrite: bool = False) -> RetroMolDuckDB:
    return RetroMolDuckDB.create(path, overwrite=overwrite)


def open_database(path: str | Path, *, read_only: bool = False) -> RetroMolDuckDB:
    return RetroMolDuckDB.open(path, read_only=read_only)
