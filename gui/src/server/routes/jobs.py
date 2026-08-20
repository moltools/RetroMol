"""Module for defining job endpoints."""

import logging
import os
import tempfile
import time
from pathlib import Path

from flask import Response, Blueprint, current_app, request, jsonify

from routes.session_store import load_session_with_items, update_item, MAX_FILE_CONTENT_BYTES
from routes.database import open_retromol_db
from routes.queue import JobStillRunningError, enqueue_and_wait
from routes.rate_limit import limiter

from retromol.chem.mol import mol_to_smiles
from retromol.chem.tagging import get_tags_mol
from retromol.model.submission import Submission
from retromol.model.rules import RuleSet
from retromol.model.readout import merge_named_paths
from retromol.model.result import Result
from retromol.pipelines.parsing import run_retromol
from retromol_fingerprint.fingerprint import TOKEN_LINK
from retromol_synthesis.reconstruction import reconstruct_linear_readout

from retromol_antismash.io import parse_antismash_gbk, AntiSmashOptions
from retromol_antismash.modules import LinearReadout, bgc_primary_sequence, linear_readout, ModuleType
from retromol_antismash.inference.registry import annotate_region
from retromol_antismash.inference.model_paras import ParasModel

from routes.discovery import get_discovery_context

blp_search_compound = Blueprint("search_compound", __name__)
blp_submit_compound = Blueprint("submit_compound", __name__)
blp_reconstruct_compound = Blueprint("reconstruct_compound", __name__)
blp_submit_gene_cluster = Blueprint("submit_gene_cluster", __name__)
blp_get_cluster_readout = Blueprint("get_cluster_readout", __name__)
blp_reconstruct_gene_cluster = Blueprint("reconstruct_gene_cluster", __name__)

# Task functions below run inside an RQ worker process, not a Flask request context --
# current_app.logger isn't available there, so they use a plain module logger instead.
# Route handlers (still running in-request) keep using current_app.logger as before.
logger = logging.getLogger(__name__)


DEFAULT_LIMIT = 10
MAX_LIMIT = 50

# PARAS' own default (retromol_antismash.inference.model_paras.ParasModel.threshold) --
# repeated here since submit_gene_cluster now builds a ParasModel explicitly (see
# _build_paras_model) rather than relying on ParasModel's constructor default.
PARAS_THRESHOLD_DEFAULT = 0.1


def _build_paras_model(threshold: float) -> ParasModel:
    """
    Build a PARAS A-domain substrate specificity model configured with the given
    probability threshold.

    A fresh instance is built per request rather than registered globally (contrast
    with `retromol_antismash.inference.registry.register_domain_model`, which keeps
    whichever instance was registered first) since the threshold is now a per-upload
    choice -- see `annotate_region`'s `domain_models` parameter. This is cheap:
    the actual model file is only downloaded/loaded lazily, the first time
    `.predict()` runs, and that underlying load is itself cached (see
    `retromol_antismash.inference.model_paras.load_paras_model`).

    :param threshold: probability threshold below which a prediction is discarded
    :return: a configured ParasModel instance
    """
    return ParasModel(
        threshold=threshold,
        model_path=os.getenv("PARAS_MODEL_PATH"),
        cache_dir=os.getenv("PARAS_CACHE_DIR", "paras_cache"),
    )


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
                    es.entry_id AS id,
                    es.name AS name,
                    es.database_name AS database_name,
                    es.url AS url,
                    e.raw AS raw
                FROM entry_sources es
                JOIN entries e ON e.id = es.entry_id
                WHERE e.type = 'compound'
                    AND e.raw IS NOT NULL
                    AND lower(es.name) LIKE lower(?)
                ORDER BY
                    CASE
                        WHEN lower(es.name) = lower(?) THEN 0
                        WHEN lower(es.name) LIKE lower(?) THEN 1
                        ELSE 2
                    END,
                    es.name
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
            "databaseName": database_name,
            "databaseIdentifier": entry_id,
            "url": url,
        }
        for entry_id, name, database_name, url, raw in rows
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


def run_compound_parse(
    session_id: str, item_id: str, name: str | None, smiles: str, match_stereochemistry: bool
) -> tuple[dict, int]:
    """
    Parse a compound with RetroMol and persist the result onto its session item.

    Runs as an RQ job (see routes/queue.py). Marks the item "processing" itself,
    once a worker actually starts this job -- not when the HTTP request was
    received -- so the SSE-driven status badge reflects real queue wait time rather
    than jumping to "processing" the instant a request lands, even if the queue is
    backed up and no worker has picked it up yet.

    :param session_id: the session the item belongs to
    :param item_id: the item to update
    :param name: display name to (re)apply to the item
    :param smiles: the SMILES to parse
    :param match_stereochemistry: whether to match stereochemistry-aware rules
    :return: a (response body, HTTP status code) pair
    """
    t0 = time.time()

    def mark_processing(it: dict) -> None:
        it["name"] = name or it.get("name")
        it["smiles"] = smiles or it.get("smiles")
        _set_item_status_inplace(it, "processing")

    if not update_item(session_id, item_id, mark_processing):
        logger.warning("run_compound_parse: failed to mark item as processing: %s", item_id)
        return {"error": "Item not found during update"}, 404

    try:
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
        logger.exception("run_compound_parse: error for item_id=%s", item_id)

        def mark_error(it: dict) -> None:
            _set_item_status_inplace(it, "error", error_message=str(e))

        update_item(session_id, item_id, mark_error)

        elapsed = int((time.time() - t0) * 1000)
        return {"ok": False, "status": "error", "elapsed_ms": elapsed, "error": str(e)}, 500

    elapsed = int((time.time() - t0) * 1000)
    logger.info("run_compound_parse: finished item_id=%s elapsed_ms=%s", item_id, elapsed)

    return {"ok": True, "status": "done", "elapsed_ms": elapsed}, 200


@blp_submit_compound.post("/api/submitCompound")
@limiter.limit("60 per minute")
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

    try:
        body, status = enqueue_and_wait(run_compound_parse, session_id, item_id, name, smiles, match_stereochemistry)
    except JobStillRunningError as e:
        return jsonify({"error": str(e)}), 503

    return jsonify(body), status


DB_MATCHING_SEQUENCE_NOTE = (
    "This is the exact primary sequence the database-population pipeline stores for "
    "this molecule (read directly off result.linear_readout with no backbone "
    "reconstruction, no eligibility/orientation filtering, no length filtering: "
    "every path found in the assembly graph is included, tailoring events like "
    "glycosylation/methylation among them). Every path is merged into this one "
    "sequence, longest first, ties broken lexicographically, joined by \"<LINK>\". "
    "The same merge database/scripts/common.py's primary_sequence_from_result "
    "applies before storing a compound. Use this to search, not a reconstructed "
    "candidate above: reconstruct_linear_readout builds a separate identified-only "
    "readout to attempt backbone reconstruction, which diverges from what's actually "
    "stored for a real fraction of compounds (anything with unidentified content). "
    "Querying with that instead silently caps how similar a match can ever score, "
    "even against the molecule's own database entry. To search with just one "
    "biosynthetic chain rather than the whole molecule, split this sequence on "
    "\"<LINK>\" and use one of the resulting subsequences instead."
)


def _db_matching_primary_sequences(result: Result) -> list[dict]:
    """
    The single primary sequence read directly off result.linear_readout -- the same
    representation database/scripts/common.py's primary_sequence_from_result uses to
    populate the persistent database (itself reproducing the original
    db_scripts/create_database.py recipe, recovered from git history at commit
    72a1bbb). Every path is merged into one sequence, nothing filtered out (see
    retromol.model.readout.LinearReadout.primary_sequence / merge_named_paths).
    Shaped identically to a Reconstruction dict, in a one-element list, so the
    frontend can render it the same way it does reconstruct_linear_readout's
    candidates, just without ever calling that function -- see
    DB_MATCHING_SEQUENCE_NOTE for why the two aren't interchangeable as a database
    query.

    :param result: a parsed RetroMol Result
    :return: a one-element list holding the single merged Reconstruction-shaped
        dict, or empty if the readout has no paths at all
    """
    tagged_input_smiles = mol_to_smiles(result.submission.mol, include_tags=True)

    named_paths = [
        [[name, list(get_tags_mol(node.mol))] for name, node in zip(
            [n.identity.matched_rule.name if n.is_identified else "X" for n in path], path
        )]
        for path in result.linear_readout.paths
    ]

    if not named_paths:
        return []

    primary_sequence = merge_named_paths(named_paths, key=lambda item: item[0], link_item=[TOKEN_LINK, []])

    return [{
        "tagged_input_smiles": tagged_input_smiles,
        "tagged_backbone_smiles": None,
        "primary_sequence": primary_sequence,
        "backbone_warning": DB_MATCHING_SEQUENCE_NOTE,
        "ordered": True,
    }]


def run_compound_reconstruction(item_payload: dict | None) -> tuple[dict, int]:
    """
    Reconstruct candidate primary sequences from an already-parsed compound Result.

    Runs as an RQ job (see routes/queue.py).

    :param item_payload: the session item's stored RetroMol Result payload
    :return: a (response body, HTTP status code) pair
    """
    try:
        result = Result.from_dict(item_payload)
        reconstructions = reconstruct_linear_readout(result)
        reconstructions_as_dicts = [rec.to_dict() for rec in reconstructions]
        db_matching_sequences = _db_matching_primary_sequences(result)
        return {
            "ok": True,
            "status": "done",
            "data": reconstructions_as_dicts,
            "dbMatchingSequences": db_matching_sequences,
        }, 200

    except Exception as e:
        logger.exception("run_compound_reconstruction: failed")
        return {"ok": False, "error": f"Failed to reconstruct compound: {e}"}, 500


@blp_reconstruct_compound.post("/api/reconstructCompound")
@limiter.limit("60 per minute")
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
        body, status = enqueue_and_wait(run_compound_reconstruction, item.get("payload"))
    except JobStillRunningError as e:
        return jsonify({"error": str(e)}), 503

    return jsonify(body), status


def run_gene_cluster_parse(
    session_id: str, item_id: str, name: str | None, file_content: str, paras_threshold: float
) -> tuple[dict, int]:
    """
    Parse an antiSMASH GenBank file into a linear module readout and persist it.

    Runs as an RQ job (see routes/queue.py); see run_compound_parse for why status
    marking happens here rather than in the route handler.

    :param session_id: the session the item belongs to
    :param item_id: the item to update
    :param name: display name to (re)apply to the item
    :param file_content: the raw antiSMASH GenBank file contents
    :param paras_threshold: minimum PARAS prediction probability to call a substrate
    :return: a (response body, HTTP status code) pair
    """
    t0 = time.time()

    def mark_processing(it: dict) -> None:
        it["name"] = name or it.get("name")
        it["fileContent"] = file_content
        it["parasThreshold"] = paras_threshold
        _set_item_status_inplace(it, "processing")

    if not update_item(session_id, item_id, mark_processing):
        logger.warning("run_gene_cluster_parse: failed to mark item as processing: %s", item_id)
        return {"error": "Item not found during update"}, 404

    try:
        paras_model = _build_paras_model(paras_threshold)

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
            # No gene-level (PFAM/HMM) model is run here, since gene-level
            # classification isn't part of this feature -- gene_models defaults to
            # whatever's globally registered, which is nothing.
            annotate_region(region, domain_models=[paras_model])
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
            it["parasThreshold"] = paras_threshold
            it["payload"] = result_payload
            if score is not None:
                it["score"] = score
            _set_item_status_inplace(it, "done")

        update_item(session_id, item_id, mark_done)

    except Exception as e:
        logger.exception("run_gene_cluster_parse: error for item_id=%s", item_id)

        def mark_error(it: dict) -> None:
            _set_item_status_inplace(it, "error", error_message=str(e))

        update_item(session_id, item_id, mark_error)

        elapsed = int((time.time() - t0) * 1000)
        return {"ok": False, "status": "error", "elapsed_ms": elapsed, "error": str(e)}, 500

    elapsed = int((time.time() - t0) * 1000)
    return {"ok": True, "status": "done", "elapsed_ms": elapsed}, 200


@blp_submit_gene_cluster.post("/api/submitGeneCluster")
@limiter.limit("60 per minute")
def submit_gene_cluster() -> tuple[Response, int]:
    """
    Endpoint to submit a gene cluster for processing.

    `parasThreshold` (optional, default PARAS_THRESHOLD_DEFAULT) sets the minimum
    prediction probability PARAS requires to call a substrate for an NRPS A-domain
    module -- lower surfaces more (lower-confidence) predictions, higher keeps only
    the most confident ones. Doesn't affect PKS modules, which are classified
    directly from domain anatomy rather than through PARAS.
    """
    payload = request.get_json(force=True) or {}

    session_id = payload.get("sessionId")
    item_id = payload.get("itemId")
    name = payload.get("name")
    file_content = payload.get("fileContent")
    paras_threshold = payload.get("parasThreshold", PARAS_THRESHOLD_DEFAULT)

    if not session_id or not item_id:
        return jsonify({"error": "Missing sessionId or itemId"}), 400

    if not isinstance(file_content, str) or not file_content:
        return jsonify({"error": "fileContent is required"}), 400

    if len(file_content.encode("utf-8")) > MAX_FILE_CONTENT_BYTES:
        max_mb = MAX_FILE_CONTENT_BYTES // (1024 * 1024)
        return jsonify({"error": f"fileContent exceeds the {max_mb} MB limit"}), 400

    if (
        not isinstance(paras_threshold, (int, float))
        or isinstance(paras_threshold, bool)
        or not (0.0 <= paras_threshold <= 1.0)
    ):
        return jsonify({"error": "parasThreshold must be a number between 0 and 1"}), 400
    paras_threshold = float(paras_threshold)

    full_sess = load_session_with_items(session_id)
    if full_sess is None:
        return jsonify({"error": "Session not found"}), 404

    item = next((it for it in full_sess.get("items", []) if it.get("id") == item_id), None)
    if item is None:
        return jsonify({"error": "Item not found"}), 404

    if item.get("kind") != "cluster":
        return jsonify({"error": "Item is not a gene cluster"}), 400

    try:
        body, status = enqueue_and_wait(run_gene_cluster_parse, session_id, item_id, name, file_content, paras_threshold)
    except JobStillRunningError as e:
        return jsonify({"error": str(e)}), 503

    return jsonify(body), status


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

        A plain Redis dict read -- no compute to isolate, so this stays direct
        rather than going through the job queue.
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


def run_gene_cluster_reconstruction(item_payload: dict | None) -> tuple[dict, int]:
    """
    Map a gene cluster's parsed linear module readout(s) onto the same
    matching-rule vocabulary a compound's primary sequence is drawn from.

    Runs as an RQ job (see routes/queue.py).

    :param item_payload: the session item's stored antiSMASH readouts payload
    :return: a (response body, HTTP status code) pair
    """
    try:
        ctx = get_discovery_context()
        readouts_data = (item_payload or {}).get("readouts") or []

        data = []
        for readout_data in readouts_data:
            readout = LinearReadout.from_dict(readout_data)
            names, _tokens = bgc_primary_sequence(readout, ctx.ruleset)
            data.append({"id": readout.id, "primary_sequence": names})

        return {"ok": True, "status": "done", "data": data}, 200

    except Exception as e:
        logger.exception("run_gene_cluster_reconstruction: failed")
        return {"ok": False, "error": str(e)}, 404


@blp_reconstruct_gene_cluster.post("/api/reconstructGeneCluster")
@limiter.limit("60 per minute")
def reconstruct_gene_cluster() -> tuple[Response, int]:
    """
    Endpoint for mapping a gene cluster's parsed linear module readout(s) onto the
    same matching-rule vocabulary a compound's primary sequence is drawn from -- the
    gene-cluster equivalent of `/api/reconstructCompound`, so a BGC upload can seed a
    discovery query the same way a compound reconstruction does.

    One entry is returned per antiSMASH region readout, each already in biosynthetic
    order (see `retromol_antismash.modules.bgc_primary_sequence`). Unlike a compound's
    reconstruction, there's no per-region edited-sequence override yet.
    """
    payload = request.get_json(force=True) or {}

    session_id = payload.get("sessionId")
    item_id = payload.get("itemId")

    full_sess = load_session_with_items(session_id)
    if full_sess is None:
        current_app.logger.warning(f"reconstruct_gene_cluster: session not found: {session_id}")
        return jsonify({"error": "Session not found"}), 404

    item = next((it for it in full_sess.get("items", []) if it.get("id") == item_id), None)
    if item is None:
        current_app.logger.warning(f"reconstruct_gene_cluster: item not found: {item_id}")
        return jsonify({"error": "Item not found"}), 404

    if item.get("kind") != "cluster":
        current_app.logger.warning(f"reconstruct_gene_cluster: wrong kind={item.get('kind')}")
        return jsonify({"error": "Item is not a gene cluster"}), 400

    try:
        body, status = enqueue_and_wait(run_gene_cluster_reconstruction, item.get("payload"))
    except JobStillRunningError as e:
        return jsonify({"error": str(e)}), 503

    return jsonify(body), status
