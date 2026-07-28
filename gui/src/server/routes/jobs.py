"""Module for defining job endpoints."""

import os
import tempfile
import time
from pathlib import Path

from flask import Response, Blueprint, current_app, request, jsonify

from routes.session_store import load_session_with_items, update_item, MAX_FILE_CONTENT_BYTES
from routes.database import open_retromol_db

from retromol.model.submission import Submission
from retromol.model.rules import RuleSet
from retromol.model.result import Result
from retromol.pipelines.parsing import run_retromol
from retromol_synthesis.reconstruction import reconstruct_linear_readout

from retromol_antismash.io import parse_antismash_gbk, AntiSmashOptions
from retromol_antismash.modules import linear_readout, ModuleType
from retromol_antismash.inference.registry import annotate_region, register_domain_model
from retromol_antismash.inference.model_paras import ParasModel

blp_search_compound = Blueprint("search_compound", __name__)
blp_submit_compound = Blueprint("submit_compound", __name__)
blp_reconstruct_compound = Blueprint("reconstruct_compound", __name__)
blp_submit_gene_cluster = Blueprint("submit_gene_cluster", __name__)
blp_get_cluster_readout = Blueprint("get_cluster_readout", __name__)


DEFAULT_LIMIT = 10
MAX_LIMIT = 50


def _ensure_paras_model_registered() -> None:
    """
    Register the PARAS A-domain substrate specificity model, if not already registered.

    .. note::
        `register_domain_model` no-ops if a model with the same name is already
        registered, so this is safe to call on every request. Constructing
        `ParasModel` here is cheap -- the actual model file is only downloaded/loaded
        lazily, the first time `.predict()` runs. No gene-level (PFAM/HMM) model is
        registered here, since gene-level classification isn't part of this feature.
    """
    register_domain_model(ParasModel(
        model_path=os.getenv("PARAS_MODEL_PATH"),
        cache_dir=os.getenv("PARAS_CACHE_DIR", "paras_cache"),
    ))


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
        _ensure_paras_model_registered()

        tmp_path: Path | None = None
        try:
            with tempfile.NamedTemporaryFile(mode="w", suffix=".gbk", delete=False) as tmp:
                tmp.write(file_content)
                tmp_path = Path(tmp.name)

            regions = parse_antismash_gbk(tmp_path, AntiSmashOptions())
        finally:
            if tmp_path is not None:
                tmp_path.unlink(missing_ok=True)

        if not regions:
            raise ValueError(
                "No antiSMASH regions found in this GenBank file. "
                "Make sure it is a GenBank file annotated by antiSMASH."
            )

        readouts = []
        for region in regions:
            annotate_region(region)
            readouts.append(linear_readout(region))

        # Rough analog of a compound's parse "coverage" score: fraction of modules
        # with a confident identification. PKS modules are classified directly from
        # domain anatomy (no ML model involved), so they always count as confident;
        # NRPS (A-domain) modules only count if PARAS resolved a non-"unknown"
        # substrate. Without this, a PKS-only cluster would score 0%/blank even
        # though every module parsed correctly, since PARAS never runs on PKS.
        all_modules = [m for r in readouts for m in r.modules]
        confident_modules = [
            m for m in all_modules
            if m.type == ModuleType.PKS
            or (m.predicted_substrate and m.predicted_substrate.name != "unknown")
        ]
        score = (len(confident_modules) / len(all_modules)) if all_modules else None

        result_payload = {"readouts": [r.to_dict() for r in readouts]}

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


@blp_get_cluster_readout.post("/api/getClusterReadout")
def get_cluster_readout() -> tuple[Response, int]:
    """
    Endpoint for fetching a gene cluster's parsed linear module readout(s).

    .. note::
        `payload` is deliberately stripped from items before they're sent to the
        client as part of the session (see `strip_property_from_dict` in
        session_store.py), the same way a compound's `payload` never reaches the
        client directly -- only via `/api/reconstructCompound`. This is the
        gene-cluster equivalent of that endpoint.
    """
    payload = request.get_json(force=True) or {}

    session_id = payload.get("sessionId")
    item_id = payload.get("itemId")

    full_sess = load_session_with_items(session_id)
    if full_sess is None:
        return jsonify({"error": "Session not found"}), 404

    item = next((it for it in full_sess.get("items", []) if it.get("id") == item_id), None)
    if item is None:
        return jsonify({"error": "Item not found"}), 404

    if item.get("kind") != "cluster":
        return jsonify({"error": "Item is not a gene cluster"}), 400

    data = item.get("payload") or {"readouts": []}
    return jsonify({"ok": True, "status": "done", "data": data}), 200
