from retromol_alignment.aligner import PairwiseAligner
from retromol_alignment.pairwise import Converter, align


def rerank(
    query: list[str],
    targets: list[list[str]],
    aligner: PairwiseAligner,
    converter: Converter,
) -> list[tuple[float, bool]]:
    scores: list[tuple[float, bool]] = []

    for target in targets:
        score1, _, _ = align(aligner, query, target, converter)
        score2, _, _ = align(aligner, query, target[::-1], converter)
        score = max(score1, score2)

        inverted = False
        if score2 > score1:
            inverted = True

        scores.append((score, inverted))

    return scores