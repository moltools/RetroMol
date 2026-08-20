"""Enrichment tab: search up to MAX_ENTRY_SEARCH_RESULTS db entries, select some of
them, and test whether the selection is enriched for any annotation term (phylogeny,
chemical class, ...) compared to its background -- entries of the same type(s) not
selected. One Fisher's exact test per term observed on the selection, Benjamini-Hochberg
corrected across all of them.
"""

from flask import Blueprint, Response, jsonify, request
from scipy.stats import fisher_exact

from retromol_database.duckdb import ENTRY_TYPES
from routes.database import open_retromol_db

blp_entry_search = Blueprint("entry_search", __name__)
blp_enrichment_analysis = Blueprint("enrichment_analysis", __name__)

# Mirrors the frontend's client-side guardrail; must also be enforced here since the
# frontend check is UX-only (see session_store.py's MAX_SESSION_ITEMS for the same pattern).
MAX_ENTRY_SEARCH_RESULTS = 100
MAX_ENRICHMENT_SELECTION = 100


def _entry_payload(entry) -> dict:
    return {
        "id": entry.id,
        "name": entry.name,
        "url": entry.url,
        "type": entry.type,
        "sources": [
            {"name": s.name, "databaseName": s.database_name, "url": s.url} for s in entry.sources
        ],
    }


@blp_entry_search.get("/api/entrySearch")
def entry_search() -> tuple[Response, int]:
    """
    Search entries by name/id, for the Enrichment tab's "query and select" step.

    :return: a tuple containing the search results and an HTTP status code
    """
    query = (request.args.get("q") or "").strip()
    if not query:
        return jsonify({"error": "q is required"}), 400

    entry_type = request.args.get("type")
    if entry_type == "all":
        entry_type = None
    if entry_type is not None and entry_type not in ENTRY_TYPES:
        return jsonify({"error": f"type must be one of {ENTRY_TYPES} or 'all'"}), 400

    limit = request.args.get("limit", default=MAX_ENTRY_SEARCH_RESULTS, type=int)
    if limit is None or not (1 <= limit <= MAX_ENTRY_SEARCH_RESULTS):
        return jsonify({"error": f"limit must be an integer between 1 and {MAX_ENTRY_SEARCH_RESULTS}"}), 400

    try:
        with open_retromol_db() as db:
            results = db.search_entries(query, entry_type=entry_type, limit=limit)
    except Exception as e:
        return jsonify({"error": str(e)}), 503

    return jsonify({"results": [_entry_payload(e) for e in results]}), 200


def _benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Benjamini-Hochberg FDR correction. Returns q-values in the same order as `p_values`."""
    m = len(p_values)
    if m == 0:
        return []

    order = sorted(range(m), key=lambda i: p_values[i])
    q_values = [0.0] * m

    running_min = 1.0
    for rank, idx in reversed(list(enumerate(order, start=1))):
        q = p_values[idx] * m / rank
        running_min = min(running_min, q)
        q_values[idx] = min(running_min, 1.0)

    return q_values


@blp_enrichment_analysis.post("/api/enrichmentAnalysis")
def enrichment_analysis() -> tuple[Response, int]:
    """
    Test whether a selected set of entries is enriched (or depleted) for any annotation
    term compared to its background -- entries of the same type(s), excluding the
    selection itself. One two-sided Fisher's exact test per term observed on at least
    one selected entry, Benjamini-Hochberg corrected across all tested terms.

    :return: a tuple containing the enrichment results and an HTTP status code
    """
    body = request.get_json(silent=True) or {}
    entry_ids = body.get("entryIds")

    if not isinstance(entry_ids, list) or not entry_ids or not all(isinstance(x, str) and x for x in entry_ids):
        return jsonify({"error": "entryIds must be a non-empty list of non-empty strings"}), 400

    entry_ids = list(dict.fromkeys(entry_ids))  # de-dupe, preserve order
    if len(entry_ids) > MAX_ENRICHMENT_SELECTION:
        return jsonify({"error": f"entryIds cannot have more than {MAX_ENRICHMENT_SELECTION} entries"}), 400

    try:
        with open_retromol_db() as db:
            selected_entries = db.get_entries(entry_ids)
            found_ids = {e.id for e in selected_entries}
            missing_ids = [i for i in entry_ids if i not in found_ids]
            if missing_ids:
                return jsonify({"error": f"unknown entry id(s): {missing_ids[:5]}"}), 400

            background_types = sorted({e.type for e in selected_entries})
            background_total = db.count_entries_by_type(background_types)
            n_selected = len(entry_ids)

            selected_counts = db.annotation_term_counts(entry_ids)
            background_counts = db.annotation_term_counts_for_types(background_types)
            terms = db.annotation_terms_by_ids(list(selected_counts.keys()))

            rows = []
            p_values = []
            for term_id, a in selected_counts.items():
                term_total = background_counts.get(term_id, a)
                b = n_selected - a
                c = term_total - a
                d = (background_total - n_selected) - c

                _, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")

                expected_rate = term_total / background_total if background_total else 0.0
                observed_rate = a / n_selected if n_selected else 0.0
                fold_enrichment = (observed_rate / expected_rate) if expected_rate > 0 else None

                term = terms.get(term_id)
                rows.append({
                    "termId": term_id,
                    "category": term.category if term else None,
                    "rank": term.rank if term else None,
                    "label": term.label if term else term_id,
                    "selectedWithTerm": a,
                    "selectedTotal": n_selected,
                    "backgroundWithTerm": term_total,
                    "backgroundTotal": background_total,
                    "foldEnrichment": fold_enrichment,
                    "direction": "enriched" if fold_enrichment is not None and fold_enrichment > 1 else "depleted",
                    "pValue": float(p_value),
                })
                p_values.append(float(p_value))

            q_values = _benjamini_hochberg(p_values)
            for row, q in zip(rows, q_values):
                row["qValue"] = q

            rows.sort(key=lambda r: r["qValue"])
    except Exception as e:
        return jsonify({"error": str(e)}), 503

    return jsonify({"results": rows}), 200
