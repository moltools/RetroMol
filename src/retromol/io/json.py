"""This module provides JSON I/O functions for the RetroMol package."""

import json
from typing import Any, Generator

import ijson


def iter_json(path: str, jsonl: bool = False) -> Generator[Any, None, None]:
    """
    Stream items from a JSON array or a JSON Lines (JSONL) file.

    :param path: path to the JSON or JSONL file
    :param jsonl: if True, treat the file as JSONL (one JSON object per line). If False, assume a single JSON array
    :yield: parsed JSON objects
    """
    with open(path, "rb") as f:
        if jsonl:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                yield json.loads(line)
        else:
            yield from ijson.items(f, "item")
