"""This module provides JSON I/O functions for the RetroMol package."""

import json
import gzip
from pathlib import Path
from typing import Any, Generator

import ijson


def open_json_file(path: str | Path, mode: str = "rt"):
    """Open plain or gzip JSON/JSONL using its .gz suffix, with bounded buffers.

    gzip level 6 balances file size and writing speed. Text files use UTF-8.
    Append mode supports concatenated gzip members, which gzip readers handle
    transparently. Callers must close the file (prefer a context manager).
    """
    if "b" not in mode and "t" not in mode:
        mode += "t"
    kwargs = {} if "b" in mode else {"encoding": "utf-8"}
    if str(path).lower().endswith(".gz"):
        return gzip.open(path, mode, compresslevel=6, **kwargs)
    return open(path, mode, **kwargs)


def dumps_compact(data: Any) -> str:
    """JSON without formatting whitespace; preserve all values and Unicode."""
    return json.dumps(data, separators=(",", ":"), ensure_ascii=False)


def iter_json(path: str | Path, jsonl: bool | None = None) -> Generator[Any, None, None]:
    """
    Stream items from a plain or gzip JSON array or JSON Lines (JSONL) file.

    :param path: Path to the JSON or JSONL file.
    :param jsonl: True for JSONL, False for a JSON array; by default infer JSONL
        from .jsonl/.ndjson (optionally followed by .gz).
    :yield: Parsed JSON objects.
    """
    if jsonl is None:
        name = str(path).lower().removesuffix(".gz")
        jsonl = name.endswith((".jsonl", ".ndjson"))
    with open_json_file(path, "rb") as f:
        if jsonl:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                yield json.loads(line)
        else:
            yield from ijson.items(f, "item")
