"""Browse tab: list/filter database entries together with their annotations, and
export the same (filtered) set as a TSV file -- compound smiles/inchikey or BGC
accession/id, alongside whatever phylogeny/chemical_class annotations are on record.
"""

import csv
import io

from flask import Blueprint, Response, jsonify, request

from retromol_database.duckdb import BrowseEntry, ENTRY_TYPES
from routes.database import open_retromol_db

blp_browse_entries = Blueprint("browse_entries", __name__)
blp_annotation_terms = Blueprint("annotation_terms", __name__)
blp_export_entries = Blueprint("export_entries", __name__)

TSV_COLUMNS = [
    "id", "type", "name", "smiles", "inchikey", "url", "sources",
    "phylogeny_type", "genus", "species", "chemical_classes",
]


def _entry_type_param() -> tuple[str | None, str | None]:
    """Read+validate the shared `type` query param. Returns (entry_type, error)."""
    entry_type = request.args.get("type")
    if entry_type in (None, "", "all"):
        return None, None
    if entry_type not in ENTRY_TYPES:
        return None, f"type must be one of {ENTRY_TYPES} or 'all'"
    return entry_type, None


def _browse_entry_payload(entry: BrowseEntry) -> dict:
    return {
        "id": entry.id,
        "type": entry.type,
        "name": entry.name,
        "url": entry.url,
        "smiles": entry.raw if entry.type == "compound" else None,
        "sources": [
            {"name": s.name, "databaseName": s.database_name, "url": s.url} for s in entry.sources
        ],
        "phylogenyType": entry.phylogeny_type,
        "genus": entry.genus,
        "species": entry.species,
        "chemicalClasses": entry.chemical_classes,
    }


def _tsv_row(entry: BrowseEntry) -> list[str]:
    return [
        entry.id,
        entry.type,
        entry.name,
        entry.raw if entry.type == "compound" else "",
        entry.id if entry.type == "compound" else "",  # entries are keyed by InChIKey for compounds
        entry.url or "",
        ";".join(f"{s.database_name}:{s.name}" for s in entry.sources),
        entry.phylogeny_type or "",
        entry.genus or "",
        entry.species or "",
        ";".join(entry.chemical_classes),
    ]


@blp_annotation_terms.get("/api/annotationTerms")
def annotation_terms() -> tuple[Response, int]:
    """
    Every annotation term in the database, for the Browse tab's filter dropdown.

    :return: a tuple containing the terms payload and an HTTP status code
    """
    category = request.args.get("category") or None

    try:
        with open_retromol_db() as db:
            terms = db.list_annotation_terms(category=category)
    except Exception as e:
        return jsonify({"error": str(e)}), 503

    return jsonify({
        "terms": [
            {"id": t.id, "category": t.category, "rank": t.rank, "label": t.label, "parentId": t.parent_id}
            for t in terms
        ]
    }), 200


@blp_browse_entries.get("/api/browseEntries")
def browse_entries() -> tuple[Response, int]:
    """
    List (optionally filtered) database entries with their annotations, for the
    Browse tab's table. Not paginated server-side -- see module docstring; the
    frontend paginates client-side over the full filtered set.

    :return: a tuple containing the entries payload and an HTTP status code
    """
    entry_type, error = _entry_type_param()
    if error:
        return jsonify({"error": error}), 400

    term_id = request.args.get("termId") or None

    try:
        with open_retromol_db() as db:
            entries = db.browse_entries(entry_type=entry_type, term_id=term_id)
    except Exception as e:
        return jsonify({"error": str(e)}), 503

    return jsonify({"entries": [_browse_entry_payload(e) for e in entries]}), 200


@blp_export_entries.get("/api/browseEntries.tsv")
def export_entries_tsv() -> Response:
    """
    Export the same (optionally filtered) set `browseEntries` returns as a TSV
    file download.

    :return: a TSV file response
    """
    entry_type, error = _entry_type_param()
    if error:
        return jsonify({"error": error}), 400

    term_id = request.args.get("termId") or None

    try:
        with open_retromol_db() as db:
            entries = db.browse_entries(entry_type=entry_type, term_id=term_id)
    except Exception as e:
        return jsonify({"error": str(e)}), 503

    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n")
    writer.writerow(TSV_COLUMNS)
    for entry in entries:
        writer.writerow(_tsv_row(entry))

    return Response(
        buf.getvalue(),
        mimetype="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="retromol_entries.tsv"'},
    )
