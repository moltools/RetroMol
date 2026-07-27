"""Module for defining job endpoints."""

import time

from flask import Response, Blueprint, current_app, request, jsonify

from routes.session_store import load_session_with_items, update_item, MAX_FILE_CONTENT_BYTES
from routes.database import open_retromol_db

from retromol.model.submission import Submission
from retromol.model.rules import RuleSet
from retromol.model.result import Result
from retromol.pipelines.parsing import run_retromol
from retromol_synthesis.reconstruction import reconstruct_linear_readout

blp_search_compound = Blueprint("search_compound", __name__)
blp_submit_compound = Blueprint("submit_compound", __name__)
blp_reconstruct_compound = Blueprint("reconstruct_compound", __name__)
blp_submit_gene_cluster = Blueprint("submit_gene_cluster", __name__)


DEFAULT_LIMIT = 10
MAX_LIMIT = 50


@blp_search_compound.get("/api/searchCompound")
def search_compound_by_name():
    """
    Autocomplete endpoint for compounds by name-like query.
    """
    q = (request.args.get("q") or "").strip()
    if not q:
        return jsonify({"rows": [], "rowCount": 0}), 200

    try:
        limit = int(request.args.get("limit", DEFAULT_LIMIT))
    except ValueError:
        limit = DEFAULT_LIMIT

    limit = max(1, min(MAX_LIMIT, limit))
    like = f"%{q}%"

    try:
        with open_retromol_db() as db:
            rows = db.con.execute(
                """
                SELECT
                    min(id) AS id,
                    name,
                    url,
                    raw
                FROM entries
                WHERE type = 'compound'
                    AND raw IS NOT NULL
                    AND lower(name) LIKE lower(?)
                GROUP BY name, url, raw
                ORDER BY
                    CASE
                        WHEN lower(name) = lower(?) THEN 0
                        WHEN lower(name) LIKE lower(?) THEN 1
                        ELSE 2
                    END,
                    name
                LIMIT ?
                """,
                [like, q, f"{q}%", limit],
            ).fetchall()

    except Exception as e:
        current_app.logger.exception("search_compound_by_name: DuckDB query failed")
        return jsonify({"error": str(e), "rows": [], "rowCount": 0}), 500

    out = [
        {
            "name": name,
            "smiles": raw,
            "databaseName": "RetroMol",
            "databaseIdentifier": entry_id,
            "url": url,
        }
        for entry_id, name, url, raw in rows
        if name and raw
    ]

    return jsonify({"rows": out, "rowCount": len(out)}), 200


def _set_item_status_inplace(item: dict, status: str, error_message: str | None = None) -> None:
    """
    Update the status and error message of an item in place.

    :param item: the item dictionary to update
    :param status: the new status string
    :param error_message: optional error message string
    """
    item["status"] = status
    item["updatedAt"] = int(time.time() * 1000)

    if error_message is not None:
        item["errorMessage"] = error_message
    else:
        if "errorMessage" in item:
            item["errorMessage"] = None


@blp_submit_compound.post("/api/submitCompound")
def submit_compound() -> tuple[Response, int]:
    """
    Endpoint to submit a compound for processing.
    """
    payload = request.get_json(force=True) or {}

    session_id = payload.get("sessionId")
    item_id = payload.get("itemId")
    name = payload.get("name")
    smiles = payload.get("smiles", None)
    match_stereochemistry = payload.get("matchStereochemistry", False)

    if smiles is None or not isinstance(smiles, str):
        current_app.logger.warning(f"submit_compound: smiles is not a string")
        return jsonify({"error": "smiles is not a string"}), 400

    current_app.logger.info(f"submit_compound called: session_id={session_id} item_id={item_id}")

    if not session_id or not item_id:
        current_app.logger.warning("submit_compound: missing sessionId or itemId")
        return jsonify({"error": "Missing sessionId or itemId"}), 400

    # Validate session + item exists and kind is correct
    full_sess = load_session_with_items(session_id)
    if full_sess is None:
        current_app.logger.warning(f"submit_compound: session not found: {session_id}")
        return jsonify({"error": "Session not found"}), 404

    item = next((it for it in full_sess.get("items", []) if it.get("id") == item_id), None)
    if item is None:
        current_app.logger.warning(f"submit_compound: item not found: {item_id}")
        return jsonify({"error": "Item not found"}), 404

    if item.get("kind") != "compound":
        current_app.logger.warning(f"submit_compound: wrong kind={item.get('kind')}")
        return jsonify({"error": "Item is not a compound"}), 400

    t0 = time.time()

    # Set status=processing early on this item only
    def mark_processing(it: dict) -> None:
        """
        Update item details and mark as processing.

        :param it: the item dictionary to update
        """
        it["name"] = name or it.get("name")
        it["smiles"] = smiles or it.get("smiles")
        _set_item_status_inplace(it, "processing")

    ok = update_item(session_id, item_id, mark_processing)
    if not ok:
        current_app.logger.warning(f"submit_compound: failed to mark item as processing: {item_id}")
        return jsonify({"error": "Item not found during update"}), 404

    try:
        # Heavy work
        submission = Submission(name="app_submission", smiles=smiles)
        rules = RuleSet.load_default(match_stereochemistry=match_stereochemistry)
        result = run_retromol(submission=submission, rules=rules)
        coverage = result.calculate_coverage()
        result_as_dict = result.to_dict()

        def mark_done(it: dict) -> None:
            it["name"] = name or it.get("name")
            it["smiles"] = smiles or it.get("smiles")
            it["matchStereochemistry"] = match_stereochemistry
            it["score"] = coverage
            it["payload"] = result_as_dict

            _set_item_status_inplace(it, "done")

        update_item(session_id, item_id, mark_done)

    except Exception as e:
        current_app.logger.exception(f"submit_compound: error for item_id={item_id}")

        def mark_error(it: dict) -> None:
            _set_item_status_inplace(it, "error", error_message=str(e))

        update_item(session_id, item_id, mark_error)

        elapsed = int((time.time() - t0) * 1000)
        return jsonify({
            "ok": False,
            "status": "error",
            "elapsed_ms": elapsed,
            "error": str(e),
        }), 500

    elapsed = int((time.time() - t0) * 1000)
    current_app.logger.info(f"submit_compound: finished item_id={item_id} elapsed_ms={elapsed}")

    return jsonify({
        "ok": True,
        "status": "done",
        "elapsed_ms": elapsed,
    }), 200


@blp_reconstruct_compound.post("/api/reconstructCompound")
def reconstruct_compound() -> tuple[Response, int]:
    """
    Endpoint for reconstructing a compound from a RetroMol result.
    """
    payload = request.get_json(force=True) or {}

    session_id = payload.get("sessionId")
    item_id = payload.get("itemId")

    full_sess = load_session_with_items(session_id)
    if full_sess is None:
        current_app.logger.warning(f"reconstruct_compound: session not found: {session_id}")
        return jsonify({"error": "Session not found"}), 404

    item = next((it for it in full_sess.get("items", []) if it.get("id") == item_id), None)
    if item is None:
        current_app.logger.warning(f"submit_compound: item not found: {item_id}")
        return jsonify({"error": "Item not found"}), 404

    if item.get("kind") != "compound":
        current_app.logger.warning(f"submit_compound: wrong kind={item.get('kind')}")
        return jsonify({"error": "Item is not a compound"}), 400

    try:
        result = Result.from_dict(item["payload"])
        reconstructions = reconstruct_linear_readout(result)
        reconstructions_as_dicts = [rec.to_dict() for rec in reconstructions]
        return jsonify({"ok": True, "status": "done", "data": reconstructions_as_dicts}), 200

    except Exception as e:
        current_app.logger.exception(f"submit_compound: error for item_id={item_id}")
        return jsonify({"ok": False, "error": "Item not found during update"}), 404


@blp_submit_gene_cluster.post("/api/submitGeneCluster")
def submit_gene_cluster() -> tuple[Response, int]:
    """
    Endpoint to submit a gene cluster for processing.
    """
    payload = request.get_json(force=True) or {}

    session_id = payload.get("sessionId")
    item_id = payload.get("itemId")
    name = payload.get("name")
    file_content = payload.get("fileContent")

    if not session_id or not item_id:
        return jsonify({"error": "Missing sessionId or itemId"}), 400

    if not isinstance(file_content, str) or not file_content:
        return jsonify({"error": "fileContent is required"}), 400

    if len(file_content.encode("utf-8")) > MAX_FILE_CONTENT_BYTES:
        max_mb = MAX_FILE_CONTENT_BYTES // (1024 * 1024)
        return jsonify({"error": f"fileContent exceeds the {max_mb} MB limit"}), 400

    full_sess = load_session_with_items(session_id)
    if full_sess is None:
        return jsonify({"error": "Session not found"}), 404

    item = next((it for it in full_sess.get("items", []) if it.get("id") == item_id), None)
    if item is None:
        return jsonify({"error": "Item not found"}), 404

    if item.get("kind") != "cluster":
        return jsonify({"error": "Item is not a gene cluster"}), 400

    t0 = time.time()

    def mark_processing(it: dict) -> None:
        it["name"] = name or it.get("name")
        it["fileContent"] = file_content
        _set_item_status_inplace(it, "processing")

    if not update_item(session_id, item_id, mark_processing):
        return jsonify({"error": "Item not found during update"}), 404

    try:
        # Fill this in with your BGC parsing / fingerprinting / scoring logic.
        result_payload = {}
        score = None

        def mark_done(it: dict) -> None:
            it["name"] = name or it.get("name")
            it["fileContent"] = file_content
            it["payload"] = result_payload
            if score is not None:
                it["score"] = score
            _set_item_status_inplace(it, "done")

        update_item(session_id, item_id, mark_done)

    except Exception as e:
        current_app.logger.exception(f"submit_gene_cluster: error for item_id={item_id}")

        def mark_error(it: dict) -> None:
            _set_item_status_inplace(it, "error", error_message=str(e))

        update_item(session_id, item_id, mark_error)

        elapsed = int((time.time() - t0) * 1000)
        return jsonify({"ok": False, "status": "error", "elapsed_ms": elapsed, "error": str(e)}), 500

    elapsed = int((time.time() - t0) * 1000)
    return jsonify({"ok": True, "status": "done", "elapsed_ms": elapsed}), 200
