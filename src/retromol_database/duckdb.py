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


@dataclass(frozen=True)
class BrowseEntry:
    id: str
    type: EntryType
    name: str
    url: str | None
    raw: str | None
    sources: list[EntrySource]
    phylogeny_type: str | None
    genus: str | None
    species: str | None
    chemical_classes: list[str]


@dataclass(frozen=True)
class AnnotationStats:
    with_annotation_count: int
    without_annotation_count: int
    counts_by_category: list[Count]
    phylogeny_type_counts: list[Count]
    top_genera: list[Count]
    chemical_class_counts: list[Count]


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
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS annotation_terms (
                id VARCHAR PRIMARY KEY,
                category VARCHAR NOT NULL,
                rank VARCHAR,
                label VARCHAR NOT NULL,
                parent_id VARCHAR
            )
            """
        )
        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS entry_annotations (
                entry_id VARCHAR NOT NULL,
                term_id VARCHAR NOT NULL,
                PRIMARY KEY (entry_id, term_id)
            )
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

    def add_annotation_term(
        self,
        *,
        term_id: str,
        category: str,
        label: str,
        rank: str | None = None,
        parent_id: str | None = None,
    ) -> str:
        """Add (or no-op if already present) a single annotation term. `term_id` is a
        caller-supplied deterministic slug (e.g. "phylogeny:genus:bacterium:streptomyces")
        so repeated calls for the same term across many entries are idempotent."""
        self.con.execute(
            """
            INSERT INTO annotation_terms (id, category, rank, label, parent_id)
            VALUES (?, ?, ?, ?, ?)
            ON CONFLICT (id) DO NOTHING
            """,
            [term_id, category, rank, label, parent_id],
        )
        return term_id

    def link_entry_annotation(self, entry_id: str, term_id: str) -> None:
        self.con.execute(
            """
            INSERT INTO entry_annotations (entry_id, term_id)
            VALUES (?, ?)
            ON CONFLICT (entry_id, term_id) DO NOTHING
            """,
            [entry_id, term_id],
        )

    def add_phylogeny_annotation(
        self,
        entry_id: str,
        *,
        type_label: str | None,
        genus: str | None,
        species: str | None,
    ) -> None:
        """Link `entry_id` to whichever of type/genus/species is resolvable, linking every
        level (not just the most specific one) so term-level queries at any rank don't need
        to walk the parent_id chain. `species` is ignored if `genus` isn't given."""
        if not type_label:
            return

        type_id = self.add_annotation_term(
            term_id=f"phylogeny:type:{type_label.lower()}",
            category="phylogeny",
            rank="type",
            label=type_label,
        )
        self.link_entry_annotation(entry_id, type_id)

        if not genus:
            return

        genus_id = self.add_annotation_term(
            term_id=f"phylogeny:genus:{type_label.lower()}:{genus.lower()}",
            category="phylogeny",
            rank="genus",
            label=genus,
            parent_id=type_id,
        )
        self.link_entry_annotation(entry_id, genus_id)

        if not species:
            return

        species_id = self.add_annotation_term(
            term_id=f"phylogeny:species:{type_label.lower()}:{genus.lower()}:{species.lower()}",
            category="phylogeny",
            rank="species",
            label=species,
            parent_id=genus_id,
        )
        self.link_entry_annotation(entry_id, species_id)

    def add_flat_annotation(self, entry_id: str, *, category: str, label: str) -> None:
        """Link `entry_id` to a single, non-hierarchical term (e.g. chemical_class, bioactivity)."""
        term_id = self.add_annotation_term(
            term_id=f"{category}:{label.lower()}",
            category=category,
            label=label,
        )
        self.link_entry_annotation(entry_id, term_id)

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
            SELECT id, category, rank, label, parent_id
            FROM annotation_terms
            WHERE id IN (SELECT UNNEST(?))
            """,
            [list(term_ids)],
        ).fetchall()
        return {
            str(row[0]): AnnotationTerm(
                id=str(row[0]), category=str(row[1]), rank=row[2], label=str(row[3]), parent_id=row[4]
            )
            for row in rows
        }

    def list_annotation_terms(self, category: str | None = None) -> list[AnnotationTerm]:
        """Every annotation term in the database (not just ones with entries counted),
        for populating a filter dropdown in the Browse tab."""
        where_sql = ""
        params: list[object] = []
        if category is not None:
            where_sql = "WHERE category = ?"
            params.append(category)

        rows = self.con.execute(
            f"""
            SELECT id, category, rank, label, parent_id
            FROM annotation_terms
            {where_sql}
            ORDER BY category, rank NULLS FIRST, label
            """,
            params,
        ).fetchall()
        return [
            AnnotationTerm(id=str(r[0]), category=str(r[1]), rank=r[2], label=str(r[3]), parent_id=r[4])
            for r in rows
        ]

    def _browse_where_clause(
        self, *, entry_type: str | None, term_id: str | None
    ) -> tuple[str, list[object]]:
        where_sql = []
        params: list[object] = []

        if entry_type is not None:
            entry_type = _normalize_entry_type(entry_type)
            where_sql.append("e.type = ?")
            params.append(entry_type)

        if term_id is not None:
            where_sql.append("e.id IN (SELECT entry_id FROM entry_annotations WHERE term_id = ?)")
            params.append(term_id)

        where_clause = f"WHERE {' AND '.join(where_sql)}" if where_sql else ""
        return where_clause, params

    def count_browse_entries(self, *, entry_type: str | None = None, term_id: str | None = None) -> int:
        """Total number of entries a `browse_entries(...)` call with the same filters would
        match, regardless of its `limit` -- lets callers warn when a result set was truncated."""
        where_clause, params = self._browse_where_clause(entry_type=entry_type, term_id=term_id)
        return int(
            self.con.execute(f"SELECT count(*) FROM entries e {where_clause}", params).fetchone()[0]
        )

    def browse_entries(
        self, *, entry_type: str | None = None, term_id: str | None = None, limit: int | None = None
    ) -> list[BrowseEntry]:
        """Every entry (optionally filtered by type and/or a single annotation term id --
        matching phylogeny at any rank or a chemical class), with its sources and
        annotations attached. Used for both the Browse tab's table and its TSV export.

        `limit` bounds how many rows are fetched -- callers should always pass one (see
        MAX_BROWSE_ENTRIES in routes/browse.py) since an unfiltered call over a
        multi-million-row database would otherwise materialize the whole table in memory.
        Use `count_browse_entries` with the same filters to tell whether the result was
        truncated. Results are ordered by id, so `limit` always returns the same prefix.
        """
        if limit is not None and limit < 1:
            raise ValueError("limit must be >= 1")

        where_clause, params = self._browse_where_clause(entry_type=entry_type, term_id=term_id)

        limit_clause = ""
        if limit is not None:
            limit_clause = "LIMIT ?"
            params = [*params, limit]

        rows = self.con.execute(
            f"""
            SELECT e.id, e.raw, e.type, e.primary_sequence, e.fingerprint
            FROM entries e
            {where_clause}
            ORDER BY e.id
            {limit_clause}
            """,
            params,
        ).fetchall()

        entry_ids = [str(row[0]) for row in rows]
        sources_by_id = self._sources_for_entry_ids(entry_ids)

        annotation_rows = (
            self.con.execute(
                """
                SELECT ea.entry_id, t.category, t.rank, t.label
                FROM entry_annotations ea
                JOIN annotation_terms t ON t.id = ea.term_id
                WHERE ea.entry_id IN (SELECT UNNEST(?))
                """,
                [entry_ids],
            ).fetchall()
            if entry_ids
            else []
        )

        phylogeny_type_by_id: dict[str, str] = {}
        genus_by_id: dict[str, str] = {}
        species_by_id: dict[str, str] = {}
        chemical_classes_by_id: dict[str, list[str]] = {}
        for entry_id, category, rank, label in annotation_rows:
            entry_id = str(entry_id)
            if category == "phylogeny" and rank == "type":
                phylogeny_type_by_id[entry_id] = str(label)
            elif category == "phylogeny" and rank == "genus":
                genus_by_id[entry_id] = str(label)
            elif category == "phylogeny" and rank == "species":
                species_by_id[entry_id] = str(label)
            elif category == "chemical_class":
                chemical_classes_by_id.setdefault(entry_id, []).append(str(label))

        out: list[BrowseEntry] = []
        for row in rows:
            entry_id = str(row[0])
            entry = _entry_from_row(row, sources_by_id.get(entry_id, []))
            out.append(
                BrowseEntry(
                    id=entry.id,
                    type=entry.type,
                    name=entry.name,
                    url=entry.url,
                    raw=entry.raw,
                    sources=entry.sources,
                    phylogeny_type=phylogeny_type_by_id.get(entry_id),
                    genus=genus_by_id.get(entry_id),
                    species=species_by_id.get(entry_id),
                    chemical_classes=chemical_classes_by_id.get(entry_id, []),
                )
            )
        return out

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

        phylogeny_type_rows = self.con.execute(
            """
            SELECT t.label, count(DISTINCT ea.entry_id)
            FROM entry_annotations ea
            JOIN annotation_terms t ON t.id = ea.term_id
            WHERE t.category = 'phylogeny' AND t.rank = 'type'
            GROUP BY t.label
            ORDER BY count(DISTINCT ea.entry_id) DESC
            """
        ).fetchall()
        phylogeny_type_counts = [Count(label=str(row[0]), count=int(row[1])) for row in phylogeny_type_rows]

        genus_rows = self.con.execute(
            """
            SELECT t.label, count(DISTINCT ea.entry_id)
            FROM entry_annotations ea
            JOIN annotation_terms t ON t.id = ea.term_id
            WHERE t.category = 'phylogeny' AND t.rank = 'genus'
            GROUP BY t.label
            ORDER BY count(DISTINCT ea.entry_id) DESC
            LIMIT 15
            """
        ).fetchall()
        top_genera = [Count(label=str(row[0]), count=int(row[1])) for row in genus_rows]

        chem_class_rows = self.con.execute(
            """
            SELECT t.label, count(DISTINCT ea.entry_id)
            FROM entry_annotations ea
            JOIN annotation_terms t ON t.id = ea.term_id
            WHERE t.category = 'chemical_class'
            GROUP BY t.label
            ORDER BY count(DISTINCT ea.entry_id) DESC
            """
        ).fetchall()
        chemical_class_counts = [Count(label=str(row[0]), count=int(row[1])) for row in chem_class_rows]

        return AnnotationStats(
            with_annotation_count=with_annotation_count,
            without_annotation_count=without_annotation_count,
            counts_by_category=counts_by_category,
            phylogeny_type_counts=phylogeny_type_counts,
            top_genera=top_genera,
            chemical_class_counts=chemical_class_counts,
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
