"""Module for defining job endpoints."""

import time

from flask import Response, Blueprint, current_app, request, jsonify

from routes.helpers import get_unique_identifier
from routes.session_store import load_session_with_items, update_item

from retromol.pipelines.parsing import run_retromol

blp_submit_compound = Blueprint("submit_compound", __name__)
blp_submit_gene_cluster = Blueprint("submit_gene_cluster", __name__)


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
    smiles = payload.get("smiles")

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
        # # Heavy work
        # generator = _setup_fingerprint_generator()
        # (
        #     tagged_smiles,
        #     coverages,
        #     fp_hex_strings,
        #     linear_readout
        # ) = _compute_compound(generator, smiles)
        #
        # # Set final status=done and store results on this item only
        # def mark_done(it: dict) -> None:
        #     it["name"] = name or it.get("name")
        #     it["smiles"] = smiles or it.get("smiles")
        #     it["taggedSmiles"] = tagged_smiles
        #     it["retrofingerprints"] = [
        #         {
        #             "id": get_unique_identifier(),
        #             "retrofingerprint512": fp_hex,
        #             "score": cov,
        #         }
        #         for cov, fp_hex in zip(coverages, fp_hex_strings, strict=True)
        #     ]
        #     it["primarySequences"] = linear_readout
        #     _set_item_status_inplace(it, "done")
        #
        # update_item(session_id, item_id, mark_done)

        # TODO: users need to be able to supply their own parsing rules, supplying own monomers should disallow fingerprinting though
        # TODO: users should be able to view the reaction graph
        # TODO: users should be able to view the 'compound->linear view->primary sequence' visualization; select from reaction graph

        raise NotImplementedError("Compound parsing not available!")

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


@blp_submit_gene_cluster.post("/api/submitGeneCluster")
def submit_gene_cluster() -> tuple[Response, int]:
    """
    Endpoint to submit a gene cluster for processing.
    """
    payload = request.get_json(force=True) or {}

    t0 = time.time()

    session_id = payload.get("sessionId")
    item_id = payload.get("itemId")
    name = payload.get("name")
    file_content = payload.get("fileContent")

    current_app.logger.info(f"submit_gene_cluster called: session_id={session_id} item_id={item_id}")

    elapsed = int((time.time() - t0) * 1000)
    current_app.logger.info(f"submit_gene_cluster: finished item_id={item_id} elapsed_ms={elapsed}")

    return jsonify({"error": "Parsing gene clusters filers is currently offline!"}), 404
