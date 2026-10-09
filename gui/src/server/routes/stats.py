"""Summary statistics over the persistent database, for the workspace home tab."""

from dataclasses import asdict

from flask import Blueprint, Response, jsonify

from routes.database import open_retromol_db

blp_database_stats = Blueprint("database_stats", __name__)


def _to_camel_case(snake: str) -> str:
    head, *tail = snake.split("_")
    return head + "".join(word.capitalize() for word in tail)


def _camelize(value):
    """Recursively camelCase every dict key in `value` -- `asdict()` on a dataclass
    with nested dataclasses (e.g. AnnotationStats.coverage: list[AnnotationCoverage])
    only produces nested plain dicts/lists, so a shallow top-level-only conversion
    misses every key below the first level (that's exactly what broke
    coverage[i].with_annotation_count -> withAnnotationCount before this)."""
    if isinstance(value, dict):
        return {_to_camel_case(k): _camelize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [_camelize(v) for v in value]
    return value


@blp_database_stats.get("/api/databaseStats")
def database_stats() -> tuple[Response, int]:
    """
    Summary statistics over the entries table (counts, sequence-length distribution,
    most common building blocks, etc.), for display on the workspace home tab.

    :return: a tuple containing the stats payload and an HTTP status code
    """
    with open_retromol_db() as db:
        stats = db.stats()

    payload = _camelize(asdict(stats))
    return jsonify(payload), 200


@blp_database_stats.get("/api/annotationStats")
def annotation_stats() -> tuple[Response, int]:
    """
    Summary statistics over the annotation_terms/entry_annotations tables (phylogeny,
    chemical class, ... coverage and per-label counts), for the workspace home tab.

    :return: a tuple containing the stats payload and an HTTP status code
    """
    with open_retromol_db() as db:
        stats = db.annotation_stats()

    payload = _camelize(asdict(stats))
    return jsonify(payload), 200
