"""Aligner module for pairwise sequence alignment."""

try:
    from Bio.Align import PairwiseAligner, substitution_matrices
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


__all__ = ["PairwiseAligner", "substitution_matrices"]


def setup_aligner(
    substitution_matrix: "substitution_matrices.Array",
    *,
    mode: str = "global",
    open_internal_gap_score: float = -0.5,
    extend_internal_gap_score: float = -0.1,
    open_end_gap_score: float = -0.2,
    extend_end_gap_score: float = -0.05,
) -> "PairwiseAligner":
    """
    Setup a PairwiseAligner with the given parameters.
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("Biopython is required for using the aligner. Please install Biopython to use this functionality!")

    if mode not in ["global", "local"]:
        raise ValueError("mode must be one of 'global' or 'local'!")

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrix
    aligner.wildcard = None

    # Internal gaps: discouraged, but allowed
    aligner.open_internal_insertion_score = open_internal_gap_score
    aligner.extend_internal_insertion_score = extend_internal_gap_score
    aligner.open_internal_deletion_score = open_internal_gap_score
    aligner.extend_internal_deletion_score = extend_internal_gap_score

    # End gaps: usually slightly less bad, especially if sequences may be truncated
    aligner.open_left_insertion_score = open_end_gap_score
    aligner.extend_left_insertion_score = extend_end_gap_score
    aligner.open_right_insertion_score = open_end_gap_score
    aligner.extend_right_insertion_score = extend_end_gap_score

    aligner.open_left_deletion_score = open_end_gap_score
    aligner.extend_left_deletion_score = extend_end_gap_score
    aligner.open_right_deletion_score = open_end_gap_score
    aligner.extend_right_deletion_score = extend_end_gap_score

    return aligner
