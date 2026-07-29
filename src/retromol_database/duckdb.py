import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Literal, Sequence

import duckdb
import numpy as np

from retromol_fingerprint.fingerprint import TOKEN_UNK

ENTRY_TYPES = ("compound", "bgc")
EntryType = Literal["compound", "bgc"]
FINGERPRINT_SIZE = 1024

# Upper bound on primary-sequence length that still gets its own histogram bucket;
# anything longer is folded into a single overflow bucket so a handful of very long
# sequences can't blow up the number of bars in a chart.
STATS_SEQUENCE_LENGTH_BUCKET_MAX = 10


@dataclass(frozen=True)
class Entry:
    id: str
    name: str
    url: str | None
    raw: str | None
    type: EntryType
    primary_sequence: list[str]
    fingerprint: list[float]


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
    sequence_length_buckets: list[Count]
    top_blocks: list[Count]
    unique_block_count: int
    fully_resolved_count: int
    has_unknown_block_count: int
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


def make_entry_id(
    *,
    name: str,
    url: str | None,
    raw: str | None,
    entry_type: str,
    primary_sequence: Sequence[str],
) -> str:
    payload = {
        "name": name,
        "url": url,
        "raw": raw,
        "type": entry_type,
        "primary_sequence": list(primary_sequence),
    }
    raw = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(raw.encode("utf-8")).hexdigest()


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
        self.con.execute(
            f"""
            CREATE TABLE IF NOT EXISTS entries (
                id VARCHAR PRIMARY KEY,
                name VARCHAR NOT NULL,
                url VARCHAR,
                raw VARCHAR,
                type VARCHAR NOT NULL CHECK (type IN ('compound', 'bgc')),
                primary_sequence VARCHAR[] NOT NULL,
                fingerprint FLOAT[{FINGERPRINT_SIZE}] NOT NULL
            )
            """
        )

    def add_entry(
        self,
        *,
        name: str,
        url: str | None,
        raw: str | None,
        entry_type: str,
        primary_sequence: Sequence[str],
        fingerprint: Sequence[float] | np.ndarray,
        entry_id: str | None = None,
    ) -> str:
        entry_type = _normalize_entry_type(entry_type)
        sequence = _normalize_primary_sequence(primary_sequence)
        fp = _normalize_fingerprint(fingerprint)

        entry_id: str = entry_id or make_entry_id(
            name=name,
            url=url,
            raw=raw,
            entry_type=entry_type,
            primary_sequence=sequence,
        )

        self.con.execute(
            """
            INSERT OR REPLACE INTO entries (
                id,
                name,
                url,
                raw,
                type,
                primary_sequence,
                fingerprint
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            [entry_id, name, url, raw, entry_type, sequence, fp],
        )

        return entry_id

    def add_entries(self, entries: Iterable[Entry]) -> int:
        rows = [
            (
                entry.id,
                entry.name,
                entry.url,
                entry.raw,
                _normalize_entry_type(entry.type),
                _normalize_primary_sequence(entry.primary_sequence),
                _normalize_fingerprint(entry.fingerprint),
            )
            for entry in entries
        ]

        if not rows:
            return 0

        self.con.executemany(
            """
            INSERT OR REPLACE INTO entries (
                id,
                name,
                url,
                raw,
                type,
                primary_sequence,
                fingerprint
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
            """,
            rows,
        )

        return len(rows)

    def count(self) -> int:
        return int(self.con.execute("SELECT count(*) FROM entries").fetchone()[0])

    def stats(self, *, top_blocks_limit: int = 12) -> DatabaseStats:
        """
        Compute summary statistics over the whole entries table, for display on a
        dashboard/overview page.

        Building blocks are the tokens in each entry's primary_sequence -- e.g. amino
        acid names, PK reduction-state groups, or tailoring events like "methylation".
        TOKEN_UNK ("<UNK>") marks a block RetroMol couldn't identify and is excluded
        from top_blocks/unique_block_count since it isn't a real building block, but is
        still reflected in fully_resolved_count/has_unknown_block_count.

        :param top_blocks_limit: how many of the most frequent building blocks to return
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

        bucket_rows = self.con.execute(
            f"""
            SELECT bucket_label, count(*) AS n
            FROM (
                SELECT
                    CASE
                        WHEN len(primary_sequence) >= {STATS_SEQUENCE_LENGTH_BUCKET_MAX}
                            THEN '{STATS_SEQUENCE_LENGTH_BUCKET_MAX}+'
                        ELSE len(primary_sequence)::VARCHAR
                    END AS bucket_label,
                    least(len(primary_sequence), {STATS_SEQUENCE_LENGTH_BUCKET_MAX}) AS bucket_order
                FROM entries
            )
            GROUP BY bucket_label, bucket_order
            ORDER BY bucket_order
            """
        ).fetchall()
        sequence_length_buckets = [Count(label=str(row[0]), count=int(row[1])) for row in bucket_rows]

        top_block_rows = self.con.execute(
            """
            SELECT token, count(*) AS n
            FROM (SELECT unnest(primary_sequence) AS token FROM entries)
            WHERE token != ?
            GROUP BY token
            ORDER BY n DESC, token
            LIMIT ?
            """,
            [TOKEN_UNK, top_blocks_limit],
        ).fetchall()
        top_blocks = [Count(label=str(row[0]), count=int(row[1])) for row in top_block_rows]

        unique_block_count = int(
            self.con.execute(
                """
                SELECT count(DISTINCT token)
                FROM (SELECT unnest(primary_sequence) AS token FROM entries)
                WHERE token != ?
                """,
                [TOKEN_UNK],
            ).fetchone()[0]
        )

        resolution_row = self.con.execute(
            """
            SELECT
                count(*) FILTER (WHERE NOT list_contains(primary_sequence, ?)),
                count(*) FILTER (WHERE list_contains(primary_sequence, ?))
            FROM entries
            """,
            [TOKEN_UNK, TOKEN_UNK],
        ).fetchone()
        fully_resolved_count = int(resolution_row[0])
        has_unknown_block_count = int(resolution_row[1])

        url_row = self.con.execute(
            """
            SELECT
                count(*) FILTER (WHERE url IS NOT NULL),
                count(*) FILTER (WHERE url IS NULL)
            FROM entries
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
            sequence_length_buckets=sequence_length_buckets,
            top_blocks=top_blocks,
            unique_block_count=unique_block_count,
            fully_resolved_count=fully_resolved_count,
            has_unknown_block_count=has_unknown_block_count,
            with_source_url_count=with_source_url_count,
            without_source_url_count=without_source_url_count,
        )

    def get_entry(self, entry_id: str) -> Entry | None:
        row = self.con.execute(
            """
            SELECT id, name, url, raw, type, primary_sequence, fingerprint
            FROM entries
            WHERE id = ?
            """,
            [entry_id],
        ).fetchone()

        if row is None:
            return None

        return _entry_from_row(row)

    def iter_entries(self) -> Iterator[Entry]:
        rows = self.con.execute(
            """
            SELECT id, name, url, raw, type, primary_sequence, fingerprint
            FROM entries
            ORDER BY id
            """
        ).fetchall()

        for row in rows:
            yield _entry_from_row(row)

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
                name,
                url,
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

        return [
            SearchResult(
                entry=_entry_from_row(row[:7]),
                similarity=float(row[7]),
            )
            for row in rows
        ]

    def export_parquet(self, path: str | Path) -> None:
        self.con.execute(
            "COPY entries TO ? (FORMAT PARQUET)",
            [str(path)],
        )

def _entry_from_row(row) -> Entry:
    return Entry(
        id=str(row[0]),
        name=str(row[1]),
        url=row[2],  # don't turn into str, might be None
        raw=row[3],
        type=_normalize_entry_type(row[4]),
        primary_sequence=list(row[5]),
        fingerprint=list(row[6]),
    )


def create_database(path: str | Path, *, overwrite: bool = False) -> RetroMolDuckDB:
    return RetroMolDuckDB.create(path, overwrite=overwrite)


def open_database(path: str | Path, *, read_only: bool = False) -> RetroMolDuckDB:
    return RetroMolDuckDB.open(path, read_only=read_only)
