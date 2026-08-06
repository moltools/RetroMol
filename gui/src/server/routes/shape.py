"""PMI shape triangle: normalized principal-moments-of-inertia analysis for a batch of molecules.

Each entry's shape is computed independently, on the heavy_compute RQ queue (see
routes/queue.py) rather than in-thread -- conformer search is real CPU work, and
this server has no other background-job mechanism. Entry count and conformer budget
are both capped to keep a single job's runtime reasonable. This intentionally uses a
much lighter conformer budget than a typical offline analysis; see
MAX_NUM_CONFS/MAX_MAX_ITERS.
"""

import logging
from typing import Any

from flask import Blueprint, jsonify, request, Response

from retromol.chem.shape import compute_pmi_shape
from routes.queue import JobStillRunningError, enqueue_and_wait
from routes.rate_limit import limiter

blp_discovery_shape = Blueprint("discovery_shape", __name__)

# run_pmi_shape_batch below runs inside an RQ worker process, not a Flask request
# context -- current_app.logger isn't available there, so it uses a plain module
# logger instead.
logger = logging.getLogger(__name__)

MAX_SHAPE_ENTRIES = 30

DEFAULT_NUM_CONFS = 25
MIN_NUM_CONFS = 5
MAX_NUM_CONFS = 60

DEFAULT_MAX_ITERS = 500
MIN_MAX_ITERS = 50
MAX_MAX_ITERS = 2000


def run_pmi_shape_batch(
    ids: list[str], smiles_list: list[str], num_confs: int, max_iters: int
) -> tuple[dict, int]:
    """
    Compute normalized PMI ratios (NPR1, NPR2) for a batch of molecules, for the
    Discovery "Shape (PMI)" view's rod/disc/sphere scatter plot.

    Runs as an RQ job (see routes/queue.py). Entries with an unparseable SMILES, an
    oversized molecule, or a failed conformer search are skipped (with a reason)
    rather than failing the whole batch -- one problem molecule shouldn't block the
    rest.

    :return: a (response body, HTTP status code) pair
    """
    results: list[dict[str, Any]] = []
    skipped: list[dict[str, str]] = []

    for entry_id, smiles in zip(ids, smiles_list):
        try:
            shape = compute_pmi_shape(smiles, num_confs=num_confs, max_iters=max_iters)
        except Exception as e:
            logger.info("run_pmi_shape_batch: skipping entry_id=%s: %s", entry_id, e)
            skipped.append({"id": entry_id, "reason": str(e)})
            continue

        results.append({"id": entry_id, "npr1": shape.npr1, "npr2": shape.npr2})

    return {"ok": True, "results": results, "skipped": skipped}, 200


@blp_discovery_shape.post("/api/discoveryShape")
@limiter.limit("10 per minute")
def discovery_shape() -> tuple[Response, int]:
    """
    Validate a PMI shape request and hand the actual conformer search off to the
    heavy_compute job queue (see routes/queue.py, run_pmi_shape_batch).

    :return: a tuple containing the per-entry PMI ratios (or an error) and an HTTP status code
    """
    payload = request.get_json(force=True) or {}

    entries = payload.get("entries")
    num_confs = payload.get("numConfs", DEFAULT_NUM_CONFS)
    max_iters = payload.get("maxIters", DEFAULT_MAX_ITERS)

    if not isinstance(entries, list) or not entries:
        return jsonify({"error": "entries must be a non-empty list"}), 400

    if len(entries) > MAX_SHAPE_ENTRIES:
        return jsonify({"error": f"entries cannot have more than {MAX_SHAPE_ENTRIES} entries"}), 400

    if not isinstance(num_confs, int) or isinstance(num_confs, bool) or not (MIN_NUM_CONFS <= num_confs <= MAX_NUM_CONFS):
        return jsonify({"error": f"numConfs must be an integer between {MIN_NUM_CONFS} and {MAX_NUM_CONFS}"}), 400

    if not isinstance(max_iters, int) or isinstance(max_iters, bool) or not (MIN_MAX_ITERS <= max_iters <= MAX_MAX_ITERS):
        return jsonify({"error": f"maxIters must be an integer between {MIN_MAX_ITERS} and {MAX_MAX_ITERS}"}), 400

    ids: list[str] = []
    smiles_list: list[str] = []
    for entry in entries:
        if (
            not isinstance(entry, dict)
            or not isinstance(entry.get("id"), str)
            or not entry.get("id")
            or not isinstance(entry.get("smiles"), str)
            or not entry.get("smiles")
        ):
            return jsonify({"error": "each entry must have a non-empty 'id' and non-empty 'smiles'"}), 400
        ids.append(entry["id"])
        smiles_list.append(entry["smiles"])

    try:
        body, status = enqueue_and_wait(run_pmi_shape_batch, ids, smiles_list, num_confs, max_iters)
    except JobStillRunningError as e:
        return jsonify({"error": str(e)}), 503

    return jsonify(body), status
