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
