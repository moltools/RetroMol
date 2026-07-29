"""Summary statistics over the persistent database, for the workspace home tab."""

from dataclasses import asdict

from flask import Blueprint, Response, jsonify

from routes.database import open_retromol_db

blp_database_stats = Blueprint("database_stats", __name__)


def _to_camel_case(snake: str) -> str:
    head, *tail = snake.split("_")
    return head + "".join(word.capitalize() for word in tail)


@blp_database_stats.get("/api/databaseStats")
def database_stats() -> tuple[Response, int]:
    """
    Summary statistics over the entries table (counts, sequence-length distribution,
    most common building blocks, etc.), for display on the workspace home tab.

    :return: a tuple containing the stats payload and an HTTP status code
    """
    with open_retromol_db() as db:
        stats = db.stats()

    payload = {_to_camel_case(k): v for k, v in asdict(stats).items()}
    return jsonify(payload), 200
