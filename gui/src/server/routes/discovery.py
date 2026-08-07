"""Discovery: fingerprint-based nearest-neighbor search plus sequence alignment against the database.

The fingerprint/alignment recipe here mirrors the (no longer in the tree) script that
originally populated the production database, `db_scripts/create_database.py` — see git
history at commit 72a1bbb. In particular the vocabulary is every matching rule's name
plus its pseudonyms, and Fingerprinter uses n_hashes=2 (not the default of 1). Deviating
from this recipe would silently make query fingerprints incomparable to the ones already
stored in the database, so it must be reproduced exactly rather than "improved".
"""

import logging
import threading
import time
import uuid
from dataclasses import dataclass
from typing import Any

import numpy as np
from flask import Blueprint, Response, current_app, jsonify, request
from rdkit import Chem

from retromol.chem.fingerprint import calculate_tanimoto_similarity, mol_to_morgan_fingerprint
from retromol.model.result import Result
from retromol.model.rules import MatchingRule, RuleSet
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.msa import calculate_msa
from retromol_alignment.pairwise import Converter, align
from retromol_alignment.ranking import rerank
from retromol_alignment.scoring import HARDCODED_PK_SCORING, create_tanimoto_scoring_matrix
from retromol_antismash.modules import LinearReadout as BgcLinearReadout, bgc_primary_sequence
from retromol_database.duckdb import FINGERPRINT_SIZE, Entry, SearchResult
from retromol_fingerprint.fingerprint import TOKEN_UNK, Fingerprinter, Vocabulary
from retromol_synthesis.reconstruction import reconstruct_linear_readout

from routes.database import open_retromol_db
from routes.queue import JobStillRunningError, enqueue_and_wait, enqueue_job
from routes.rate_limit import limiter
from routes.session_store import load_session_with_items, save_item, update_item
from routes.shape import DEFAULT_MAX_ITERS, DEFAULT_NUM_CONFS, MAX_SHAPE_ENTRIES, run_pmi_shape_batch

blp_discovery_monomer_names = Blueprint("discovery_monomer_names", __name__)
blp_discovery_query = Blueprint("discovery_query", __name__)
blp_discovery_msa = Blueprint("discovery_msa", __name__)
blp_discovery_compare_tanimoto = Blueprint("discovery_compare_tanimoto", __name__)
blp_submit_discovery_query = Blueprint("submit_discovery_query", __name__)
blp_get_discovery_query_result = Blueprint("get_discovery_query_result", __name__)

# Task functions below (run_discovery_query, run_discovery_msa, run_tanimoto_batch)
# run inside an RQ worker process, not a Flask request context -- current_app.logger
# isn't available there, so they use a plain module logger instead. Route handlers
# (still running in-request) keep using current_app.logger as before.
logger = logging.getLogger(__name__)


MAX_N = 1000
MAX_TOP_X = 200
# Query plus up to MAX_TOP_X candidates, with a little headroom -- this endpoint is
# meant to be called with sequences already returned by /api/discoveryQuery, which
# already caps candidate count at MAX_TOP_X.
MAX_MSA_SEQUENCES = 250
ENTRY_TYPES_FOR_QUERY = ("compound", "bgc", "both")

# Same headroom as MAX_MSA_SEQUENCES -- this endpoint is also meant to be called with
# candidates already returned by /api/discoveryQuery (capped at MAX_TOP_X).
MAX_COMPARE_ENTRIES = 250
MIN_MORGAN_RADIUS = 1
MAX_MORGAN_RADIUS = 6
MIN_MORGAN_NBITS = 16
MAX_MORGAN_NBITS = 8192

# Radius/fingerprint size used when precomputing the Compare metric at query-submit
# time (see run_discovery_query_job) -- mirrors discovery_compare_tanimoto's own
# defaults (payload.get("radius", 2) / payload.get("nBits", 2048)) so the precomputed
# value matches what a user would get from the interactive picker's default state.
DEFAULT_COMPARE_RADIUS = 2
DEFAULT_COMPARE_NBITS = 2048

# A saved discovery-query session item -- see run_discovery_query_job/submit_discovery_query.
DISCOVERY_QUERY_KIND = "discoveryQuery"
MAX_DISCOVERY_QUERY_ITEMS = 5

# What a candidate's alignment score gets normalized against, as a percentage.
# "subsequence" only considers the query's own length, so a short query that
# matches well as a subsequence of a much longer candidate can still score near
# 100% -- intentional for that mode, since it's meant to surface subsequence
# matches. "longest_sequence" instead normalizes against whichever of the two
# sequences is longer, which penalizes exactly that case and favors matches
# between similarly-sized sequences. To add a new mode: add its name here and a
# branch in discovery_query() below.
ALIGNMENT_SCORE_MODES = ("subsequence", "longest_sequence")



# Group-level PKS pseudo-tokens: a BGC's PKS module only ever resolves to a reduction
# level (see retromol_antismash.modules.PKSExtenderUnit), never a specific side chain,
# so it can't be named as one specific matching rule the way an NRPS module can. These
# aren't matching-rule names themselves -- they're the pseudonym every rule at that
# reduction level shares (already part of the fingerprint vocabulary via `pseudonyms`;
# see `_per_monomer_tokens`) -- but they must also be added to the *alignment* alphabet
# for a coarse BGC call to be alignable against specific compound-side rule names at
# all. HARDCODED_PK_SCORING (retromol_alignment.scoring) pre-declares a 1.0 similarity
# between each group token and every rule name at that reduction level; every other
# pairing (including between different group tokens) defaults to 0 similarity.
PK_GROUP_TOKENS = ("PK_A", "PK_B", "PK_C", "PK_D")


@dataclass(frozen=True)
class DiscoveryContext:
    """Everything expensive to (re)build: the ruleset, fingerprint vocab, and alignment scoring matrix."""

    ruleset: RuleSet
    name_to_rule: dict[str, MatchingRule]
    rule_names_sorted: list[str]
    fingerprinter: Fingerprinter
    aligner: Any
    converter: Converter
    alphabet: frozenset[str]


_context: DiscoveryContext | None = None
_context_lock = threading.Lock()


def _build_context() -> DiscoveryContext:
    rules = RuleSet.load_default()

    vocab_tokens: set[str] = set()
    name_to_rule: dict[str, MatchingRule] = {}
    records: list[tuple[str, str]] = []

    for rule in rules.matching_rules:
        vocab_tokens.add(rule.name)
        vocab_tokens.update(rule.pseudonyms)
        name_to_rule.setdefault(rule.name, rule)
        records.append((rule.name, rule.smiles))

    vocab = Vocabulary(vocab_tokens)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)

    scoring_matrix = create_tanimoto_scoring_matrix(
        records,
        radius=2,
        num_bits=2048,
        stereochemistry=False,
        self_score_tokens=[TOKEN_UNK, *PK_GROUP_TOKENS],
        self_score=1.0,
        hardcoded_scores=HARDCODED_PK_SCORING,
    )
    aligner = setup_aligner(scoring_matrix, mode="global")
    alphabet = scoring_matrix.alphabet
    converter = Converter(
        to_identifier=lambda s: alphabet.index(s),
        from_identifier=lambda i: alphabet[i],
    )

    return DiscoveryContext(
        ruleset=rules,
        name_to_rule=name_to_rule,
        # PK_GROUP_TOKENS aren't matching-rule names (name_to_rule.get would miss
        # them), but they are valid primary-sequence block names -- both
        # _per_monomer_tokens and _normalize_for_alignment special-case them -- so
        # they belong in the autocomplete list the Sequence Editor searches.
        rule_names_sorted=sorted({*name_to_rule, *PK_GROUP_TOKENS}),
        fingerprinter=fingerprinter,
        aligner=aligner,
        converter=converter,
        alphabet=frozenset(alphabet),
    )


def get_discovery_context() -> DiscoveryContext:
    """
    Lazily build (once per process) and return the shared discovery context.

    Building this from scratch reparses the ruleset YAML, rebuilds every rule's RDKit
    Mol, and computes an all-pairs Tanimoto scoring matrix over every named monomer —
    expensive, so it must never happen per-request.

    :return: the shared DiscoveryContext
    """
    global _context
    if _context is not None:
        return _context
    with _context_lock:
        if _context is None:
            _context = _build_context()
    return _context


def _normalize_for_alignment(name: str, ctx: DiscoveryContext) -> str:
    """
    Map a display-layer block name to the token used in the scoring-matrix alphabet.

    "X" is how the app's existing reconstruction endpoint displays an unidentified
    block; TOKEN_UNK ("<UNK>") is what the database actually stores for the same
    concept. Any other name we don't recognize is treated the same way. Membership is
    checked against the full alphabet rather than just `name_to_rule`, since the
    alphabet also contains the PK_GROUP_TOKENS pseudo-names (e.g. "PK_A") a BGC's PKS
    module resolves to -- those aren't matching-rule names, so they'd otherwise be
    mistaken for unrecognized names and collapsed to TOKEN_UNK.

    :param name: the display-layer block name
    :param ctx: the discovery context
    :return: a token guaranteed to exist in the scoring-matrix alphabet
    """
    if name == "X" or name not in ctx.alphabet:
        return TOKEN_UNK
    return name


def _denormalize_for_display(name: str) -> str:
    """
    Map the internal unidentified-token placeholder back to the app's "X" display convention.

    :param name: an internal token (possibly TOKEN_UNK)
    :return: the display-layer name
    """
    return "X" if name == TOKEN_UNK else name


def _pair_similarity(ctx: DiscoveryContext, a: str | None, b: str | None) -> float | None:
    """
    Look up the Tanimoto similarity between two normalized (pre-display) tokens from the
    same scoring matrix used to compute the alignment -- so this is exactly "how similar
    is this substitution", not a re-derived or approximated number.

    :param ctx: the discovery context
    :param a: a normalized token, or None for a gap
    :param b: a normalized token, or None for a gap
    :return: the Tanimoto similarity in [0, 1], or None if either side is a gap
    """
    if a is None or b is None:
        return None
    return float(ctx.aligner.substitution_matrix[a, b])


def _self_alignment_score(sequence: list[str], ctx: DiscoveryContext) -> float:
    """
    Align a sequence against itself -- its own best-possible alignment score.

    :param sequence: an alignment-normalized token sequence (every token in ctx.alphabet)
    :param ctx: the discovery context
    :return: the self-alignment score
    """
    score, _, _ = align(ctx.aligner, sequence, sequence, ctx.converter)
    return score


def _normalized_pct(align_score: float, denom: float) -> float:
    """
    Convert a raw alignment score into a percentage of the best achievable score.

    :param align_score: the raw alignment score
    :param denom: the best-achievable ("self") alignment score to normalize against
    :return: a percentage clamped to [0, 100]
    """
    if denom <= 0:
        return 0.0
    return max(0.0, min(1.0, align_score / denom)) * 100.0


UPLOAD_ENTRY_ID_PREFIX = "upload:"


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    """
    Cosine similarity between two fingerprint vectors, matching the metric
    `RetroMolDuckDB.closest` uses (`array_cosine_similarity`) so upload candidates
    are ranked on the same scale as database candidates.

    :param a: first fingerprint vector
    :param b: second fingerprint vector
    :return: cosine similarity, or 0.0 if either vector has zero norm
    """
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    if denom <= 0:
        return 0.0
    return float(np.dot(a, b) / denom)


# Which uploaded item kinds count for a given requested entryType -- mirrors
# ENTRY_TYPES_FOR_QUERY, but in terms of session item "kind" rather than stored
# entry "type" (an uploaded gene cluster item has kind "cluster", stored/queried as
# entry type "bgc").
_UPLOAD_KINDS_FOR_ENTRY_TYPE: dict[str, frozenset[str]] = {
    "compound": frozenset({"compound"}),
    "bgc": frozenset({"cluster"}),
    "both": frozenset({"compound", "cluster"}),
}


def _build_compound_upload_candidate(item: dict, ctx: DiscoveryContext, query_fp: np.ndarray) -> list[SearchResult]:
    """
    Build synthetic search candidates from one uploaded (and successfully parsed)
    compound -- each of its reconstructed primary sequences becomes its own candidate,
    using whatever sequence is currently effective for it (the user's edited override
    if one was saved, otherwise the algorithm's own parse).

    :param item: the session item (kind "compound", status "done", with a payload)
    :param ctx: the discovery context
    :param query_fp: the query's fingerprint, for scoring candidates against
    :return: synthetic SearchResult candidates, one per reconstructed sequence
    """
    try:
        reconstructions = reconstruct_linear_readout(Result.from_dict(item["payload"]))
    except Exception:
        current_app.logger.exception("discovery_query: failed to reconstruct upload item_id=%s", item.get("id"))
        return []

    overrides = item.get("editedPrimarySequences") or {}
    label = item.get("name") or "Uploaded compound"

    candidates: list[SearchResult] = []
    for idx, reconstruction in enumerate(reconstructions):
        override = overrides.get(str(idx))
        effective_sequence = override if override is not None else reconstruction.to_dict()["primary_sequence"]
        names = [name for name, _tags in effective_sequence]
        if not names:
            continue

        fp = ctx.fingerprinter.encode([_per_monomer_tokens(name, ctx) for name in names])
        candidates.append(
            SearchResult(
                entry=Entry(
                    id=f"{UPLOAD_ENTRY_ID_PREFIX}{item['id']}:{idx}",
                    name=label if len(reconstructions) == 1 else f"{label} #{idx + 1}",
                    url=None,
                    # The uploaded compound's own SMILES -- same molecule for every
                    # reconstruction candidate derived from it. Carried through as `raw`
                    # so it lines up with how database compound entries store their
                    # SMILES (see the module docstring / _build_context), letting the
                    # Tanimoto compare endpoint treat both origins identically.
                    raw=item.get("smiles"),
                    type="compound",
                    primary_sequence=names,
                    fingerprint=fp.tolist(),
                ),
                similarity=_cosine_similarity(query_fp, fp),
            )
        )

    return candidates


def _build_bgc_upload_candidate(item: dict, ctx: DiscoveryContext, query_fp: np.ndarray) -> list[SearchResult]:
    """
    Build synthetic search candidates from one uploaded (and successfully parsed) gene
    cluster -- each antiSMASH region's linear readout becomes its own candidate.

    Each module's predicted substrate is mapped onto the same matching-rule vocabulary
    compound primary sequences use, via `bgc_primary_sequence` -- PKS modules resolve
    to their reduction-level pseudonym (e.g. "PK_A"), NRPS modules resolve to whichever
    rule(s) their PARAS-predicted substrate matches by structure (see
    `retromol_antismash.modules.module_primary_sequence_tokens` for why the two module
    types resolve differently). This is what makes a BGC candidate's fingerprint and
    primary sequence comparable to a compound's in the first place -- unlike compounds,
    there's no per-item edited-sequence override for gene clusters yet.

    :param item: the session item (kind "cluster", status "done", with a payload)
    :param ctx: the discovery context
    :param query_fp: the query's fingerprint, for scoring candidates against
    :return: synthetic SearchResult candidates, one per region's linear readout
    """
    readouts_data = (item.get("payload") or {}).get("readouts") or []
    label = item.get("name") or "Uploaded gene cluster"

    candidates: list[SearchResult] = []
    for idx, readout_data in enumerate(readouts_data):
        try:
            readout = BgcLinearReadout.from_dict(readout_data)
            names, tokens = bgc_primary_sequence(readout, ctx.ruleset)
        except Exception:
            current_app.logger.exception(
                "discovery_query: failed to build primary sequence for upload item_id=%s region=%s",
                item.get("id"), readout_data.get("id"),
            )
            continue

        if not names:
            continue

        fp = ctx.fingerprinter.encode(tokens)
        candidates.append(
            SearchResult(
                entry=Entry(
                    id=f"{UPLOAD_ENTRY_ID_PREFIX}{item['id']}:{idx}",
                    name=label if len(readouts_data) == 1 else f"{label} #{idx + 1}",
                    url=None,
                    raw=None,
                    type="bgc",
                    primary_sequence=names,
                    fingerprint=fp.tolist(),
                ),
                similarity=_cosine_similarity(query_fp, fp),
            )
        )

    return candidates


def _build_upload_candidates(
    session_id: str, ctx: DiscoveryContext, query_fp: np.ndarray, entry_type: str
) -> list[SearchResult]:
    """
    Build synthetic search candidates from the session's own uploaded (and successfully
    parsed) compounds and/or gene clusters, so they can compete for a spot among the
    nearest-neighbor candidates alongside the persistent database entries.

    :param session_id: the session whose uploads should be considered
    :param ctx: the discovery context
    :param query_fp: the query's fingerprint, for scoring candidates against
    :param entry_type: one of ENTRY_TYPES_FOR_QUERY, restricting which upload kinds are considered
    :return: synthetic SearchResult candidates, one per reconstructed/readout sequence
    """
    session = load_session_with_items(session_id)
    if session is None:
        current_app.logger.warning("discovery_query: sessionId not found for includeUserUploads: %s", session_id)
        return []

    allowed_kinds = _UPLOAD_KINDS_FOR_ENTRY_TYPE[entry_type]

    candidates: list[SearchResult] = []
    for item in session.get("items", []):
        kind = item.get("kind")
        if kind not in allowed_kinds or item.get("status") != "done" or not item.get("payload"):
            continue

        if kind == "compound":
            candidates.extend(_build_compound_upload_candidate(item, ctx, query_fp))
        else:
            candidates.extend(_build_bgc_upload_candidate(item, ctx, query_fp))

    return candidates


def _per_monomer_tokens(name: str, ctx: DiscoveryContext) -> list[str]:
    """
    Build the fingerprinting token list for one primary-sequence block.

    :param name: the display-layer block name
    :param ctx: the discovery context
    :return: tokens for Fingerprinter.encode (empty list -> falls back to the unknown token)
    """
    if name in PK_GROUP_TOKENS:
        # Not a matching-rule name -- the group-level pseudonym a BGC's PKS module
        # resolves to (see PK_GROUP_TOKENS). Still resolvable to fingerprint tokens
        # directly since it, and "PK", are already pseudonyms of every rule at that
        # reduction level and so already part of the fingerprint vocabulary.
        return [name, "PK"]

    rule = ctx.name_to_rule.get(name)
    if rule is None:
        return []
    tokens = {rule.name}
    tokens.update(rule.pseudonyms)
    return list(tokens)


@blp_discovery_monomer_names.get("/api/discoveryMonomerNames")
def discovery_monomer_names() -> tuple[Response, int]:
    """
    Autocomplete endpoint for valid primary-sequence block names -- matching rule
    names plus the PK_GROUP_TOKENS pseudo-names (e.g. "PK_A") used for a
    reduction-level-only PKS call.

    :return: a tuple containing a dictionary with matching names and an HTTP status code
    """
    q = (request.args.get("q") or "").strip()
    if not q:
        return jsonify({"rows": [], "rowCount": 0}), 200

    try:
        limit = int(request.args.get("limit", 10))
    except ValueError:
        limit = 10
    limit = max(1, min(50, limit))

    ctx = get_discovery_context()
    q_lower = q.lower()

    exact: list[str] = []
    prefix: list[str] = []
    contains: list[str] = []
    for name in ctx.rule_names_sorted:
        name_lower = name.lower()
        if name_lower == q_lower:
            exact.append(name)
        elif name_lower.startswith(q_lower):
            prefix.append(name)
        elif q_lower in name_lower:
            contains.append(name)

    ordered = (exact + prefix + contains)[:limit]
    return jsonify({"rows": [{"name": n} for n in ordered], "rowCount": len(ordered)}), 200


def run_discovery_query(
    primary_sequence: list[str],
    entry_type: str,
    n: int,
    top_x: int,
    score_mode: str,
    only_user_uploads: bool,
    include_user_uploads: bool,
    session_id: str | None,
) -> tuple[dict, int]:
    """
    Fingerprint a primary sequence, retrieve nearest neighbors from the database, and
    align the query against each candidate.

    Runs as an RQ job (see routes/queue.py) -- builds its own DiscoveryContext
    lazily via get_discovery_context() the first time it runs in a given worker
    process, same as any other caller of that function.

    :return: a (response body, HTTP status code) pair
    """
    ctx = get_discovery_context()

    # Fingerprint is built from the display sequence directly (Fingerprinter.encode
    # already treats an empty token list as "unknown", matching how unidentified
    # blocks were encoded when the database was originally built).
    per_monomer_tokens = [_per_monomer_tokens(name, ctx) for name in primary_sequence]
    query_fp = ctx.fingerprinter.encode(per_monomer_tokens)

    # Alignment needs every token to exist in the scoring-matrix alphabet, so
    # unidentified/unrecognized blocks collapse to the single TOKEN_UNK placeholder.
    alignment_query = [_normalize_for_alignment(name, ctx) for name in primary_sequence]

    # onlyUserUploads skips the persistent database entirely, so a user who just
    # wants to align their own uploads against each other never pays for a DB round-trip.
    db_candidates: list[SearchResult] = []
    if not only_user_uploads:
        with open_retromol_db() as db:
            if entry_type == "both":
                compound_hits = db.closest(query_fp, limit=n, entry_type="compound")
                bgc_hits = db.closest(query_fp, limit=n, entry_type="bgc")
                db_candidates = [*compound_hits, *bgc_hits]
            else:
                db_candidates = db.closest(query_fp, limit=n, entry_type=entry_type)

    # User uploads compete for a spot among the nearest N neighbors on the same
    # fingerprint-similarity footing as database entries.
    upload_candidates = (
        _build_upload_candidates(session_id, ctx, query_fp, entry_type)
        if include_user_uploads
        else []
    )

    candidates: list[SearchResult] = sorted(
        [*db_candidates, *upload_candidates], key=lambda r: r.similarity, reverse=True
    )[:n]

    if not candidates:
        return {
            "ok": True,
            "scoreMode": score_mode,
            "querySequence": [_denormalize_for_display(t) for t in alignment_query],
            "selfScore": 0.0,
            "candidatesConsidered": 0,
            "candidatesSkipped": 0,
            "results": [],
        }, 200

    try:
        self_score, _, _ = align(ctx.aligner, alignment_query, alignment_query, ctx.converter)
    except Exception as e:
        logger.exception("run_discovery_query: self-alignment failed")
        return {"error": f"Could not align the query sequence against itself: {e}"}, 400

    # Skip (rather than crash on) any candidate whose stored sequence contains a token
    # outside the current scoring-matrix alphabet -- e.g. rule data drifted since the
    # database was populated.
    usable_candidates: list[SearchResult] = []
    usable_targets: list[list[str]] = []
    skipped = 0

    for candidate in candidates:
        target = candidate.entry.primary_sequence
        if target and all(tok in ctx.alphabet for tok in target):
            usable_candidates.append(candidate)
            usable_targets.append(target)
        else:
            skipped += 1

    reranked = rerank(alignment_query, usable_targets, ctx.aligner, ctx.converter) if usable_targets else []

    # Only "longest_sequence" needs each candidate's own self-alignment score --
    # "subsequence" normalizes every candidate by the same constant (self_score),
    # which doesn't change relative ranking, so skip the extra work for it.
    target_self_scores: list[float | None]
    if score_mode == "longest_sequence":
        target_self_scores = [_self_alignment_score(target, ctx) for target in usable_targets]
    else:
        target_self_scores = [None] * len(usable_targets)

    def rank_key(align_score: float, target_self_score: float | None) -> float:
        if score_mode == "longest_sequence" and target_self_score is not None:
            return _normalized_pct(align_score, max(self_score, target_self_score))
        return align_score

    scored = list(zip(usable_candidates, usable_targets, reranked, target_self_scores))
    scored.sort(key=lambda item: rank_key(item[2][0], item[3]), reverse=True)
    top = scored[:top_x]

    results: list[dict[str, Any]] = []
    for candidate, target, (score, inverted), target_self_score in top:
        oriented_target = target[::-1] if inverted else target

        try:
            align_score, aligned_query_display, aligned_target_display = align(
                ctx.aligner, alignment_query, oriented_target, ctx.converter
            )
        except Exception:
            logger.exception(
                "run_discovery_query: alignment traceback failed for entry_id=%s", candidate.entry.id
            )
            continue

        # Self-alignment score is invariant under reversing both operands, so the
        # score computed against the un-reversed target is still valid here.
        denom = max(self_score, target_self_score) if target_self_score is not None else self_score
        normalized_pct = _normalized_pct(align_score, denom)

        aligned_similarity = [
            _pair_similarity(ctx, q_tok, t_tok) for q_tok, t_tok in zip(aligned_query_display, aligned_target_display)
        ]

        results.append({
            "entryId": candidate.entry.id,
            "name": candidate.entry.name,
            "url": candidate.entry.url,
            "type": candidate.entry.type,
            "origin": "upload" if candidate.entry.id.startswith(UPLOAD_ENTRY_ID_PREFIX) else "database",
            # The candidate's own molecular SMILES, when known -- only ever set for
            # compound entries (database compounds store it as `raw`; uploaded
            # compounds carry their own `smiles` through the same field, see
            # _build_compound_upload_candidate). None for BGCs and for any compound
            # entry that predates SMILES being stored. Lets the frontend compute a
            # real structure-based Tanimoto comparison without a second round trip
            # to look up what SMILES a given entryId corresponds to.
            "smiles": candidate.entry.raw if candidate.entry.type == "compound" else None,
            "fingerprintSimilarity": candidate.similarity,
            "primarySequence": [_denormalize_for_display(t) for t in oriented_target],
            "inverted": inverted,
            "alignmentScore": align_score,
            "normalizedAlignmentScorePct": normalized_pct,
            "alignedQuery": [
                (_denormalize_for_display(t) if t is not None else None) for t in aligned_query_display
            ],
            "alignedTarget": [
                (_denormalize_for_display(t) if t is not None else None) for t in aligned_target_display
            ],
            # Per-column Tanimoto similarity between alignedQuery[i] and alignedTarget[i],
            # None wherever either side is a gap -- lets the frontend shade each unit by
            # how good its own match was, rather than by the row's overall score.
            "alignedSimilarity": aligned_similarity,
        })

    return {
        "ok": True,
        "scoreMode": score_mode,
        "querySequence": [_denormalize_for_display(t) for t in alignment_query],
        "selfScore": self_score,
        "candidatesConsidered": len(candidates),
        "candidatesSkipped": skipped,
        "results": results,
    }, 200


@blp_discovery_query.post("/api/discoveryQuery")
@limiter.limit("60 per minute")
def discovery_query() -> tuple[Response, int]:
    """
    Validate a discovery query request and hand the actual search/alignment off to
    the heavy_compute job queue (see routes/queue.py, run_discovery_query).

    :return: a tuple containing the ranked results (or an error) and an HTTP status code
    """
    payload = request.get_json(force=True) or {}

    primary_sequence = payload.get("primarySequence")
    entry_type = payload.get("entryType")
    n = payload.get("n")
    top_x = payload.get("topX")
    score_mode = payload.get("scoreMode", "subsequence")
    only_user_uploads = bool(payload.get("onlyUserUploads", False))
    include_user_uploads = bool(payload.get("includeUserUploads", False)) or only_user_uploads
    session_id = payload.get("sessionId")

    if (
        not isinstance(primary_sequence, list)
        or not primary_sequence
        or not all(isinstance(x, str) and x for x in primary_sequence)
    ):
        return jsonify({"error": "primarySequence must be a non-empty list of non-empty strings"}), 400

    if entry_type not in ENTRY_TYPES_FOR_QUERY:
        return jsonify({"error": f"entryType must be one of {ENTRY_TYPES_FOR_QUERY}"}), 400

    if score_mode not in ALIGNMENT_SCORE_MODES:
        return jsonify({"error": f"scoreMode must be one of {ALIGNMENT_SCORE_MODES}"}), 400

    if not isinstance(n, int) or isinstance(n, bool) or not (1 <= n <= MAX_N):
        return jsonify({"error": f"n must be an integer between 1 and {MAX_N}"}), 400

    max_top_x = min(n, MAX_TOP_X)
    if not isinstance(top_x, int) or isinstance(top_x, bool) or not (1 <= top_x <= max_top_x):
        return jsonify({"error": f"topX must be an integer between 1 and {max_top_x}"}), 400

    if include_user_uploads and not (isinstance(session_id, str) and session_id):
        return jsonify({"error": "sessionId is required when includeUserUploads or onlyUserUploads is true"}), 400

    try:
        body, status = enqueue_and_wait(
            run_discovery_query,
            primary_sequence, entry_type, n, top_x, score_mode, only_user_uploads, include_user_uploads, session_id,
        )
    except JobStillRunningError as e:
        return jsonify({"error": str(e)}), 503

    return jsonify(body), status


def run_discovery_msa(
    query_sequence: list[str], ids: list[str], raw_sequences: list[list[str]]
) -> tuple[dict, int]:
    """
    Compute a multiple sequence alignment of a query against a set of sequences,
    anchored on the query as the star center (so every other sequence is ordered and
    oriented relative to it, matching how the pairwise results are already anchored).

    Runs as an RQ job (see routes/queue.py). Intended to be called with the query
    and result sequences already returned by /api/discoveryQuery -- scores for
    display are deliberately NOT recomputed here; the frontend re-attaches each
    row's existing normalizedAlignmentScorePct by id, so the number shown next to a
    row in the MSA always matches the pairwise view.

    :return: a (response body, HTTP status code) pair
    """
    ctx = get_discovery_context()

    all_ids = ["query", *ids]
    to_align = [
        [_normalize_for_alignment(name, ctx) for name in seq]
        for seq in [query_sequence, *raw_sequences]
    ]

    try:
        aligned, new_order = calculate_msa(
            ctx.aligner,
            to_align,
            ctx.converter,
            center_star=0,
            # Every sequence's orientation was already decided relative to the query
            # by the pairwise alignment step; re-flipping here would risk showing a
            # different orientation than the pairwise view already displayed.
            can_reverse=[False] * len(to_align),
        )
    except Exception as e:
        logger.exception("run_discovery_msa: MSA computation failed")
        return {"error": f"Could not compute the multiple sequence alignment: {e}"}, 400

    # center_star=0 above pins the query (to_align[0]) as the star center, so aligned[0]
    # is always the query's own row, already in the shared MSA column coordinate system --
    # every other row's per-column similarity is computed against this exact sequence.
    query_aligned_tokens = aligned[0][2]

    rows = []
    for (_, _, aligned_seq), original_index in zip(aligned, new_order):
        entry_id = all_ids[original_index]
        similarity_to_query = (
            [None] * len(aligned_seq)
            if entry_id == "query"
            else [_pair_similarity(ctx, q_tok, tok) for q_tok, tok in zip(query_aligned_tokens, aligned_seq)]
        )
        rows.append({
            "id": entry_id,
            "alignedSequence": [
                (_denormalize_for_display(t) if t is not None else None) for t in aligned_seq
            ],
            # Per-column Tanimoto similarity against the query's own row at that same MSA
            # column; all None for the query's row itself (nothing to compare it against).
            "similarityToQuery": similarity_to_query,
        })

    return {"ok": True, "rows": rows}, 200


@blp_discovery_msa.post("/api/discoveryMsa")
@limiter.limit("60 per minute")
def discovery_msa() -> tuple[Response, int]:
    """
    Validate a multiple sequence alignment request and hand the actual alignment off
    to the heavy_compute job queue (see routes/queue.py, run_discovery_msa).

    :return: a tuple containing the aligned rows (or an error) and an HTTP status code
    """
    payload = request.get_json(force=True) or {}

    query_sequence = payload.get("querySequence")
    sequences = payload.get("sequences")

    if (
        not isinstance(query_sequence, list)
        or not query_sequence
        or not all(isinstance(x, str) and x for x in query_sequence)
    ):
        return jsonify({"error": "querySequence must be a non-empty list of non-empty strings"}), 400

    if not isinstance(sequences, list) or not sequences:
        return jsonify({"error": "sequences must be a non-empty list"}), 400

    if len(sequences) > MAX_MSA_SEQUENCES:
        return jsonify({"error": f"sequences cannot have more than {MAX_MSA_SEQUENCES} entries"}), 400

    ids: list[str] = []
    raw_sequences: list[list[str]] = []
    for item in sequences:
        if (
            not isinstance(item, dict)
            or not isinstance(item.get("id"), str)
            or not item.get("id")
            or not isinstance(item.get("sequence"), list)
            or not item.get("sequence")
            or not all(isinstance(x, str) and x for x in item["sequence"])
        ):
            return jsonify({
                "error": "each entry in sequences must have a non-empty 'id' and non-empty 'sequence'"
            }), 400
        ids.append(item["id"])
        raw_sequences.append(item["sequence"])

    try:
        body, status = enqueue_and_wait(run_discovery_msa, query_sequence, ids, raw_sequences)
    except JobStillRunningError as e:
        return jsonify({"error": str(e)}), 503

    return jsonify(body), status


def run_tanimoto_batch(
    query_smiles: str, radius: int, n_bits: int, ids: list[str], smiles_list: list[str]
) -> tuple[dict, int]:
    """
    Compute a real (whole-molecule) Morgan/ECFP Tanimoto similarity between a query
    SMILES and a batch of candidate SMILES.

    Runs as an RQ job (see routes/queue.py). Unlike fingerprintSimilarity and
    normalizedAlignmentScorePct (both already returned by /api/discoveryQuery,
    computed over the coarse per-monomer-token representation used for database
    retrieval), this operates on each side's actual molecular structure.

    :return: a (response body, HTTP status code) pair
    """
    query_mol = Chem.MolFromSmiles(query_smiles)
    if query_mol is None:
        return {"error": "querySmiles could not be parsed as a valid molecule"}, 400
    query_fp = mol_to_morgan_fingerprint(query_mol, radius=radius, num_bits=n_bits)

    values: list[dict[str, Any]] = []
    invalid_count = 0
    for entry_id, smiles in zip(ids, smiles_list):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            invalid_count += 1
            continue
        fp = mol_to_morgan_fingerprint(mol, radius=radius, num_bits=n_bits)
        values.append({"id": entry_id, "tanimotoSimilarity": calculate_tanimoto_similarity(query_fp, fp)})

    return {
        "ok": True,
        "radius": radius,
        "nBits": n_bits,
        "values": values,
        "invalidCount": invalid_count,
    }, 200


@blp_discovery_compare_tanimoto.post("/api/discoveryCompareTanimoto")
@limiter.limit("10 per minute")
def discovery_compare_tanimoto() -> tuple[Response, int]:
    """
    Validate a Tanimoto comparison request and hand the actual RDKit fingerprinting
    off to the heavy_compute job queue (see routes/queue.py, run_tanimoto_batch).

    Meant for the Discovery "Compare" view, where the user can pick a radius/
    fingerprint size and see the metric recomputed on demand. Entries are meant to
    be exactly the `smiles` field already present on results from
    /api/discoveryQuery (null/omitted for BGCs and any compound missing one), pre-
    filtered by the frontend before calling this endpoint.

    :return: a tuple containing per-entry Tanimoto similarities (or an error) and an HTTP status code
    """
    payload = request.get_json(force=True) or {}

    query_smiles = payload.get("querySmiles")
    radius = payload.get("radius", 2)
    n_bits = payload.get("nBits", 2048)
    entries = payload.get("entries")

    if not isinstance(query_smiles, str) or not query_smiles.strip():
        return jsonify({"error": "querySmiles must be a non-empty string"}), 400

    if not isinstance(radius, int) or isinstance(radius, bool) or not (MIN_MORGAN_RADIUS <= radius <= MAX_MORGAN_RADIUS):
        return jsonify({"error": f"radius must be an integer between {MIN_MORGAN_RADIUS} and {MAX_MORGAN_RADIUS}"}), 400

    if not isinstance(n_bits, int) or isinstance(n_bits, bool) or not (MIN_MORGAN_NBITS <= n_bits <= MAX_MORGAN_NBITS):
        return jsonify({"error": f"nBits must be an integer between {MIN_MORGAN_NBITS} and {MAX_MORGAN_NBITS}"}), 400

    if not isinstance(entries, list) or not entries:
        return jsonify({"error": "entries must be a non-empty list"}), 400

    if len(entries) > MAX_COMPARE_ENTRIES:
        return jsonify({"error": f"entries cannot have more than {MAX_COMPARE_ENTRIES} entries"}), 400

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
        body, status = enqueue_and_wait(run_tanimoto_batch, query_smiles, radius, n_bits, ids, smiles_list)
    except JobStillRunningError as e:
        return jsonify({"error": str(e)}), 503

    return jsonify(body), status


def _set_item_status_inplace(item: dict, status: str, error_message: str | None = None) -> None:
    """
    Update the status and error message of an item in place.

    Duplicated from routes/jobs.py's identical helper rather than imported --
    jobs.py imports get_discovery_context from this module, so importing back from
    jobs.py here would create a circular import.

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


def run_discovery_query_job(session_id: str, item_id: str, settings: dict, flags: dict) -> tuple[dict, int]:
    """
    Run a full discovery query -- the base ranking plus whichever optional metrics
    were flagged at submission -- and persist the result onto its session item.

    Runs as an RQ job (see routes/queue.py), on the protected heavy_compute_pmi
    queue when flags["computePmi"] is set, otherwise the regular heavy_compute queue
    (see submit_discovery_query). Marks the item "processing" itself, once a worker
    actually starts this job, same as run_compound_parse in jobs.py.

    All requested metrics land in one `payload`, written via a single update_item
    call once everything is done -- one terminal status per query, no partial/
    incremental reveal (a PMI-flagged query stays "processing" until PMI is done too).

    :param session_id: the session the item belongs to
    :param item_id: the item to update
    :param settings: the query's own settings (primarySequence, entryType, n, topX, ...)
    :param flags: which optional metrics to precompute (computeMsa/computeCompare/computePmi)
    :return: a (response body, HTTP status code) pair
    """
    t0 = time.time()

    def mark_processing(it: dict) -> None:
        _set_item_status_inplace(it, "processing")

    if not update_item(session_id, item_id, mark_processing):
        logger.warning("run_discovery_query_job: failed to mark item as processing: %s", item_id)
        return {"error": "Item not found during update"}, 404

    try:
        query_body, query_status = run_discovery_query(
            settings["primarySequence"],
            settings["entryType"],
            settings["n"],
            settings["topX"],
            settings["scoreMode"],
            settings.get("onlyUserUploads", False),
            settings.get("includeUserUploads", False),
            session_id,
        )
        if query_status != 200:
            raise RuntimeError(query_body.get("error", "Discovery query failed"))

        results = query_body.get("results", [])
        payload: dict[str, Any] = dict(query_body)
        query_origin_smiles = settings.get("queryOriginSmiles")

        if flags.get("computeMsa") and results:
            msa_body, msa_status = run_discovery_msa(
                query_body["querySequence"],
                [r["entryId"] for r in results],
                [r["primarySequence"] for r in results],
            )
            if msa_status == 200:
                payload["msaRows"] = msa_body.get("rows", [])
            else:
                logger.warning(
                    "run_discovery_query_job: MSA precompute failed for item_id=%s: %s", item_id, msa_body
                )

        if flags.get("computeCompare") and query_origin_smiles and results:
            compare_entries = [(r["entryId"], r["smiles"]) for r in results if r.get("type") == "compound" and r.get("smiles")]
            if compare_entries:
                compare_body, compare_status = run_tanimoto_batch(
                    query_origin_smiles,
                    DEFAULT_COMPARE_RADIUS,
                    DEFAULT_COMPARE_NBITS,
                    [e[0] for e in compare_entries],
                    [e[1] for e in compare_entries],
                )
                if compare_status == 200:
                    payload["compareValues"] = compare_body.get("values", [])
                    payload["compareRadius"] = DEFAULT_COMPARE_RADIUS
                    payload["compareNBits"] = DEFAULT_COMPARE_NBITS
                else:
                    logger.warning(
                        "run_discovery_query_job: compare precompute failed for item_id=%s: %s", item_id, compare_body
                    )

        if flags.get("computePmi") and results:
            shape_entries = [(r["entryId"], r["smiles"]) for r in results if r.get("type") == "compound" and r.get("smiles")]
            if query_origin_smiles:
                shape_entries = [("query", query_origin_smiles), *shape_entries]
            shape_entries = shape_entries[:MAX_SHAPE_ENTRIES]

            if shape_entries:
                shape_body, shape_status = run_pmi_shape_batch(
                    [e[0] for e in shape_entries],
                    [e[1] for e in shape_entries],
                    DEFAULT_NUM_CONFS,
                    DEFAULT_MAX_ITERS,
                )
                if shape_status == 200:
                    payload["pmiResults"] = shape_body.get("results", [])
                    payload["pmiSkipped"] = shape_body.get("skipped", [])
                else:
                    logger.warning(
                        "run_discovery_query_job: PMI precompute failed for item_id=%s: %s", item_id, shape_body
                    )

        def mark_done(it: dict) -> None:
            it["payload"] = payload
            _set_item_status_inplace(it, "done")

        update_item(session_id, item_id, mark_done)

    except Exception as e:
        logger.exception("run_discovery_query_job: error for item_id=%s", item_id)

        def mark_error(it: dict) -> None:
            _set_item_status_inplace(it, "error", error_message=str(e))

        update_item(session_id, item_id, mark_error)

        elapsed = int((time.time() - t0) * 1000)
        return {"ok": False, "status": "error", "elapsed_ms": elapsed, "error": str(e)}, 500

    elapsed = int((time.time() - t0) * 1000)
    logger.info("run_discovery_query_job: finished item_id=%s elapsed_ms=%s", item_id, elapsed)
    return {"ok": True, "status": "done", "elapsed_ms": elapsed}, 200


@blp_submit_discovery_query.post("/api/submitDiscoveryQuery")
@limiter.limit("30 per minute")
def submit_discovery_query() -> tuple[Response, int]:
    """
    Validate a discovery query submission, persist it as a session item, and hand the
    actual computation off to an RQ worker -- fire-and-forget, unlike /api/discoveryQuery.

    Unlike the older, synchronous /api/discoveryQuery (still used internally by
    run_discovery_query_job above), this endpoint never blocks on the result: the
    created item is returned immediately with status "queued", and its status/payload
    updates flow purely through the existing SSE (item_updated) + getSession
    mechanism -- the same way compound/cluster items already work (see
    routes/jobs.py's submit_compound). A session is capped at MAX_DISCOVERY_QUERY_ITEMS
    saved queries; submitting past that is rejected rather than silently evicting one.

    :return: a tuple containing the created item and an HTTP status code
    """
    payload = request.get_json(force=True) or {}

    session_id = payload.get("sessionId")
    primary_sequence = payload.get("primarySequence")
    entry_type = payload.get("entryType")
    n = payload.get("n")
    top_x = payload.get("topX")
    score_mode = payload.get("scoreMode", "subsequence")
    only_user_uploads = bool(payload.get("onlyUserUploads", False))
    include_user_uploads = bool(payload.get("includeUserUploads", False)) or only_user_uploads
    query_origin_smiles = payload.get("queryOriginSmiles")
    flags = payload.get("flags") or {}
    name = payload.get("name") or "Discovery query"

    if not isinstance(session_id, str) or not session_id:
        return jsonify({"error": "Missing sessionId"}), 400

    if (
        not isinstance(primary_sequence, list)
        or not primary_sequence
        or not all(isinstance(x, str) and x for x in primary_sequence)
    ):
        return jsonify({"error": "primarySequence must be a non-empty list of non-empty strings"}), 400

    if entry_type not in ENTRY_TYPES_FOR_QUERY:
        return jsonify({"error": f"entryType must be one of {ENTRY_TYPES_FOR_QUERY}"}), 400

    if score_mode not in ALIGNMENT_SCORE_MODES:
        return jsonify({"error": f"scoreMode must be one of {ALIGNMENT_SCORE_MODES}"}), 400

    if not isinstance(n, int) or isinstance(n, bool) or not (1 <= n <= MAX_N):
        return jsonify({"error": f"n must be an integer between 1 and {MAX_N}"}), 400

    max_top_x = min(n, MAX_TOP_X)
    if not isinstance(top_x, int) or isinstance(top_x, bool) or not (1 <= top_x <= max_top_x):
        return jsonify({"error": f"topX must be an integer between 1 and {max_top_x}"}), 400

    if query_origin_smiles is not None and not isinstance(query_origin_smiles, str):
        return jsonify({"error": "queryOriginSmiles must be a string or null"}), 400

    if not isinstance(flags, dict):
        return jsonify({"error": "flags must be an object"}), 400

    flags_normalized = {
        "computeMsa": bool(flags.get("computeMsa", False)),
        "computeCompare": bool(flags.get("computeCompare", False)),
        "computePmi": bool(flags.get("computePmi", False)),
    }

    full_sess = load_session_with_items(session_id)
    if full_sess is None:
        return jsonify({"error": "Session not found"}), 404

    existing_query_items = [it for it in full_sess.get("items", []) if it.get("kind") == DISCOVERY_QUERY_KIND]
    if len(existing_query_items) >= MAX_DISCOVERY_QUERY_ITEMS:
        return jsonify({
            "error": f"You've reached the {MAX_DISCOVERY_QUERY_ITEMS} saved query limit -- delete one to add another"
        }), 400

    settings = {
        "primarySequence": primary_sequence,
        "entryType": entry_type,
        "scoreMode": score_mode,
        "n": n,
        "topX": top_x,
        "includeUserUploads": include_user_uploads,
        "onlyUserUploads": only_user_uploads,
        "queryOriginSmiles": query_origin_smiles,
    }

    item = {
        "id": str(uuid.uuid4()),
        "kind": DISCOVERY_QUERY_KIND,
        "name": name,
        "status": "queued",
        "errorMessage": None,
        "updatedAt": int(time.time() * 1000),
        "settings": settings,
        "flags": flags_normalized,
        "payload": None,
    }
    save_item(session_id, item)

    enqueue_job(run_discovery_query_job, session_id, item["id"], settings, flags_normalized, heavy=flags_normalized["computePmi"])

    return jsonify({"ok": True, "item": item}), 202


@blp_get_discovery_query_result.post("/api/getDiscoveryQueryResult")
def get_discovery_query_result() -> tuple[Response, int]:
    """
    Fetch a discoveryQuery session item's computed payload.

    `payload` is deliberately stripped from items before they're sent to the client
    as part of the session (see strip_property_from_dict in session_store.py, and
    routes/session.py's get_session) -- same as a compound's (only reachable via
    /api/reconstructCompound) or a cluster's (via /api/getClusterReadout). This is
    the discoveryQuery equivalent.

    A plain Redis dict read -- no compute to isolate, so this stays direct rather
    than going through the job queue.

    :return: a tuple containing the item's status/payload and an HTTP status code
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

    if item.get("kind") != DISCOVERY_QUERY_KIND:
        return jsonify({"error": "Item is not a discovery query"}), 400

    return jsonify({"ok": True, "status": item.get("status"), "payload": item.get("payload")}), 200
