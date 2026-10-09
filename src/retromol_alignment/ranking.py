import numpy as np
from scipy.optimize import linear_sum_assignment

from retromol.model.readout import split_named_paths
from retromol_fingerprint.fingerprint import TOKEN_LINK
from retromol_alignment.aligner import PairwiseAligner
from retromol_alignment.pairwise import Converter, align


def reorder_target_chains(
    query: list[str],
    target: list[str],
    aligner: PairwiseAligner,
    converter: Converter,
    link_token: str = TOKEN_LINK,
) -> tuple[float, list[str], bool]:
    """
    Reorder and reorient `target`'s chains (split on `link_token`) to best match
    `query`'s chains, order-agnostically -- a merged primary sequence's chain
    order is an artifact of how `merge_named_paths` sorted its candidate paths
    (longest first, then lexicographically), not a biosynthetic correspondence
    between two different compounds' chains, so a straight flat alignment (or
    the old whole-sequence-only `rerank`) can easily misalign chain N of one
    sequence against a structurally unrelated chain M of the other just because
    they happened to land at the same position.

    `query`'s own chain order is left untouched -- it's the fixed reference every
    target gets reordered against, the same role it already played in the old
    whole-sequence rerank (query fixed, target flipped).

    Matching is a 1:1 optimal assignment (Hungarian algorithm) over every
    (query_chain, target_chain) pair's own best-of-forward/reverse alignment
    score -- exact rather than greedy, which matters here since a greedy pick
    of the single best-scoring pair can lock out a better overall assignment
    (e.g. two chains that are each other's second-best match, but each other's
    globally best match is a chain the other also wants). Chain counts are
    small (typically 2-5), so this optimal solve is cheap.

    When the two sides have different chain counts (e.g. a tailoring event
    present on one side but not the other), the leftover chains on the larger
    side are appended after the matched ones, in their original relative order
    -- unmatched on purpose, so they end up aligning against a gap run in the
    final flat alignment (penalizing the missing content) instead of silently
    disappearing.

    :param query: alignment-normalized query token sequence (may contain link_token)
    :param target: alignment-normalized target token sequence (may contain link_token)
    :param aligner: the PairwiseAligner to score/align chain pairs with
    :param converter: Converter for the aligner's alphabet
    :param link_token: the token marking a join point between chains (default: TOKEN_LINK)
    :return: (assignment_score, reordered_target, any_chain_reversed) --
        assignment_score is the sum of each matched pair's own alignment score,
        informational/for ranking only: the score actually shown/used for a
        candidate should come from re-aligning `query` against `reordered_target`
        as one flat sequence (see routes/discovery.py's run_discovery_query), so
        chain-boundary gaps are scored consistently with every other candidate.
        any_chain_reversed is True if any matched chain needed flipping.
    """
    identity = lambda item: item  # noqa: E731 -- items here are already the plain string tokens

    query_chains = split_named_paths(query, identity, link_token)
    target_chains = split_named_paths(target, identity, link_token)

    if not query_chains or not target_chains:
        return 0.0, list(target), False

    n_q, n_t = len(query_chains), len(target_chains)
    scores = np.zeros((n_q, n_t))
    reversed_flags = np.zeros((n_q, n_t), dtype=bool)

    for i, q_chain in enumerate(query_chains):
        for j, t_chain in enumerate(target_chains):
            score_fwd, _, _ = align(aligner, q_chain, t_chain, converter)
            score_rev, _, _ = align(aligner, q_chain, t_chain[::-1], converter)
            if score_rev > score_fwd:
                scores[i, j] = score_rev
                reversed_flags[i, j] = True
            else:
                scores[i, j] = score_fwd

    # linear_sum_assignment minimizes cost; negate to maximize alignment score.
    # It handles a non-square matrix natively, matching min(n_q, n_t) pairs and
    # leaving the rest of the larger side unmatched -- exactly the "leftover
    # chains" case described above.
    row_ind, col_ind = linear_sum_assignment(-scores)

    # zip(row_ind, col_ind) is already in ascending query-chain order (that's
    # how linear_sum_assignment returns it), so the matched chains come out in
    # the query's own chain order -- the natural choice, since query is the
    # fixed reference.
    reordered: list[str] = []
    assignment_score = 0.0
    any_reversed = False

    for i, j in zip(row_ind, col_ind):
        chain = target_chains[j]
        if reversed_flags[i, j]:
            chain = chain[::-1]
            any_reversed = True

        if reordered:
            reordered.append(link_token)
        reordered.extend(chain)
        assignment_score += scores[i, j]

    matched_target_idx = set(col_ind.tolist())
    for j in range(n_t):
        if j in matched_target_idx:
            continue
        if reordered:
            reordered.append(link_token)
        reordered.extend(target_chains[j])

    return assignment_score, reordered, any_reversed
