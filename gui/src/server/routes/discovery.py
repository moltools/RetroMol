"""Discovery: fingerprint-based nearest-neighbor search plus sequence alignment against the database.

The fingerprint/alignment recipe here mirrors the (no longer in the tree) script that
originally populated the production database, `db_scripts/create_database.py` — see git
history at commit 72a1bbb. In particular the vocabulary is every matching rule's name
plus its pseudonyms, and Fingerprinter uses n_hashes=2 (not the default of 1). Deviating
from this recipe would silently make query fingerprints incomparable to the ones already
stored in the database, so it must be reproduced exactly rather than "improved".
"""

import threading
from dataclasses import dataclass
from typing import Any

from flask import Blueprint, Response, current_app, jsonify, request

from retromol.model.rules import MatchingRule, RuleSet
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.msa import calculate_msa
from retromol_alignment.pairwise import Converter, align
from retromol_alignment.ranking import rerank
from retromol_alignment.scoring import create_tanimoto_scoring_matrix
from retromol_database.duckdb import FINGERPRINT_SIZE, SearchResult
from retromol_fingerprint.fingerprint import TOKEN_UNK, Fingerprinter, Vocabulary

from routes.database import open_retromol_db

blp_discovery_monomer_names = Blueprint("discovery_monomer_names", __name__)
blp_discovery_query = Blueprint("discovery_query", __name__)
blp_discovery_msa = Blueprint("discovery_msa", __name__)


MAX_N = 1000
MAX_TOP_X = 200
# Query plus up to MAX_TOP_X candidates, with a little headroom -- this endpoint is
# meant to be called with sequences already returned by /api/discoveryQuery, which
# already caps candidate count at MAX_TOP_X.
MAX_MSA_SEQUENCES = 250
ENTRY_TYPES_FOR_QUERY = ("compound", "bgc", "both")

# What a candidate's alignment score gets normalized against, as a percentage.
# "subsequence" only considers the query's own length, so a short query that
# matches well as a subsequence of a much longer candidate can still score near
# 100% -- intentional for that mode, since it's meant to surface subsequence
# matches. "longest_sequence" instead normalizes against whichever of the two
# sequences is longer, which penalizes exactly that case and favors matches
# between similarly-sized sequences. To add a new mode: add its name here and a
# branch in discovery_query() below.
ALIGNMENT_SCORE_MODES = ("subsequence", "longest_sequence")


@dataclass(frozen=True)
class DiscoveryContext:
    """Everything expensive to (re)build: the ruleset, fingerprint vocab, and alignment scoring matrix."""

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
        self_score_tokens=[TOKEN_UNK],
        self_score=1.0,
    )
    aligner = setup_aligner(scoring_matrix, mode="global")
    alphabet = scoring_matrix.alphabet
    converter = Converter(
        to_identifier=lambda s: alphabet.index(s),
        from_identifier=lambda i: alphabet[i],
    )

    return DiscoveryContext(
        name_to_rule=name_to_rule,
        rule_names_sorted=sorted(name_to_rule),
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
    concept. Any other name we don't recognize is treated the same way.

    :param name: the display-layer block name
    :param ctx: the discovery context
    :return: a token guaranteed to exist in the scoring-matrix alphabet
    """
    if name == "X" or name not in ctx.name_to_rule:
        return TOKEN_UNK
    return name


def _denormalize_for_display(name: str) -> str:
    """
    Map the internal unidentified-token placeholder back to the app's "X" display convention.

    :param name: an internal token (possibly TOKEN_UNK)
    :return: the display-layer name
    """
    return "X" if name == TOKEN_UNK else name


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


def _per_monomer_tokens(name: str, ctx: DiscoveryContext) -> list[str]:
    """
    Build the fingerprinting token list for one primary-sequence block.

    :param name: the display-layer block name
    :param ctx: the discovery context
    :return: tokens for Fingerprinter.encode (empty list -> falls back to the unknown token)
    """
    rule = ctx.name_to_rule.get(name)
    if rule is None:
        return []
    tokens = {rule.name}
    tokens.update(rule.pseudonyms)
    return list(tokens)


@blp_discovery_monomer_names.get("/api/discoveryMonomerNames")
def discovery_monomer_names() -> tuple[Response, int]:
    """
    Autocomplete endpoint for valid primary-sequence block names (real matching rule names).

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


@blp_discovery_query.post("/api/discoveryQuery")
def discovery_query() -> tuple[Response, int]:
    """
    Fingerprint a primary sequence, retrieve nearest neighbors from the database, and
    align the query against each candidate.

    :return: a tuple containing the ranked results (or an error) and an HTTP status code
    """
    payload = request.get_json(force=True) or {}

    primary_sequence = payload.get("primarySequence")
    entry_type = payload.get("entryType")
    n = payload.get("n")
    top_x = payload.get("topX")
    score_mode = payload.get("scoreMode", "subsequence")

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

    ctx = get_discovery_context()

    # Fingerprint is built from the display sequence directly (Fingerprinter.encode
    # already treats an empty token list as "unknown", matching how unidentified
    # blocks were encoded when the database was originally built).
    per_monomer_tokens = [_per_monomer_tokens(name, ctx) for name in primary_sequence]
    query_fp = ctx.fingerprinter.encode(per_monomer_tokens)

    # Alignment needs every token to exist in the scoring-matrix alphabet, so
    # unidentified/unrecognized blocks collapse to the single TOKEN_UNK placeholder.
    alignment_query = [_normalize_for_alignment(name, ctx) for name in primary_sequence]

    with open_retromol_db() as db:
        if entry_type == "both":
            compound_hits = db.closest(query_fp, limit=n, entry_type="compound")
            bgc_hits = db.closest(query_fp, limit=n, entry_type="bgc")
            candidates: list[SearchResult] = sorted(
                [*compound_hits, *bgc_hits], key=lambda r: r.similarity, reverse=True
            )[:n]
        else:
            candidates = db.closest(query_fp, limit=n, entry_type=entry_type)

    if not candidates:
        return jsonify({
            "ok": True,
            "scoreMode": score_mode,
            "querySequence": [_denormalize_for_display(t) for t in alignment_query],
            "selfScore": 0.0,
            "candidatesConsidered": 0,
            "candidatesSkipped": 0,
            "results": [],
        }), 200

    try:
        self_score, _, _ = align(ctx.aligner, alignment_query, alignment_query, ctx.converter)
    except Exception as e:
        current_app.logger.exception("discovery_query: self-alignment failed")
        return jsonify({"error": f"Could not align the query sequence against itself: {e}"}), 400

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
            current_app.logger.exception(
                "discovery_query: alignment traceback failed for entry_id=%s", candidate.entry.id
            )
            continue

        # Self-alignment score is invariant under reversing both operands, so the
        # score computed against the un-reversed target is still valid here.
        denom = max(self_score, target_self_score) if target_self_score is not None else self_score
        normalized_pct = _normalized_pct(align_score, denom)

        results.append({
            "entryId": candidate.entry.id,
            "name": candidate.entry.name,
            "url": candidate.entry.url,
            "type": candidate.entry.type,
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
        })

    return jsonify({
        "ok": True,
        "scoreMode": score_mode,
        "querySequence": [_denormalize_for_display(t) for t in alignment_query],
        "selfScore": self_score,
        "candidatesConsidered": len(candidates),
        "candidatesSkipped": skipped,
        "results": results,
    }), 200


@blp_discovery_msa.post("/api/discoveryMsa")
def discovery_msa() -> tuple[Response, int]:
    """
    Compute a multiple sequence alignment of a query against a set of sequences,
    anchored on the query as the star center (so every other sequence is ordered and
    oriented relative to it, matching how the pairwise results are already anchored).

    Intended to be called with the query and result sequences already returned by
    /api/discoveryQuery -- scores for display are deliberately NOT recomputed here;
    the frontend re-attaches each row's existing normalizedAlignmentScorePct by id,
    so the number shown next to a row in the MSA always matches the pairwise view.

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
        current_app.logger.exception("discovery_msa: MSA computation failed")
        return jsonify({"error": f"Could not compute the multiple sequence alignment: {e}"}), 400

    rows = [
        {
            "id": all_ids[original_index],
            "alignedSequence": [
                (_denormalize_for_display(t) if t is not None else None) for t in aligned_seq
            ],
        }
        for (_, _, aligned_seq), original_index in zip(aligned, new_order)
    ]

    return jsonify({"ok": True, "rows": rows}), 200
