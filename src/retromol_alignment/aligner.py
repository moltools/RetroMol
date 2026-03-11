"""Aligner module for pairwise sequence alignment."""

try:
    from Bio.Align import PairwiseAligner, substitution_matrices
    BIOPYTHON_AVAILABLE = True
except ImportError:
    BIOPYTHON_AVAILABLE = False


__all__ = ["PairwiseAligner", "substitution_matrices"]


def setup_aligner(
    substitution_matrix: substitution_matrices.Array,
    *,
    mode: str = "global",
    open_internal_insertion_score: float = -1.0,
    extend_internal_insertion_score: float = -1.0,
    open_left_insertion_score: float = -1.0,
    extend_left_insertion_score: float = -1.0,
    open_right_insertion_score: float = -1.0,
    extend_right_insertion_score: float = -1.0,
    open_internal_deletion_score: float = -1.0,
    extend_internal_deletion_score: float = -1.0,
    open_left_deletion_score: float = -1.0,
    extend_left_deletion_score: float = -1.0,
    open_right_deletion_score: float = -1.0,
    extend_right_deletion_score: float = -1.0,
) -> PairwiseAligner:
    """
    Setup a PairwiseAligner with the given parameters.

    :param substitution_matrix: Substitution matrix to use for scoring matches/mismatches.
    :param mode: Alignment mode, either "global" or "local" (default: "global").
    :param open_internal_insertion_score: Score for opening an internal insertion (default: -1.0).
    :param extend_internal_insertion_score: Score for extending an internal insertion (default: -1.0).
    :param open_left_insertion_score: Score for opening a left insertion (default: -1.0).
    :param extend_left_insertion_score: Score for extending a left insertion (default: -1.0).
    :param open_right_insertion_score: Score for opening a right insertion (default: -1.0).
    :param extend_right_insertion_score: Score for extending a right insertion (default: -1.0).
    :param open_internal_deletion_score: Score for opening an internal deletion (default: -1.0).
    :param extend_internal_deletion_score: Score for extending an internal deletion (default: -1.0).
    :param open_left_deletion_score: Score for opening a left deletion (default: -1.0).
    :param extend_left_deletion_score: Score for extending a left deletion (default: -1.0).
    :param open_right_deletion_score: Score for opening a right deletion (default: -1.0).
    :param extend_right_deletion_score: Score for extending a right deletion (default: -1.0).
    :return: PairwiseAligner object.
    :raises ValueError: If mode is not "global" or "local".
    :raises ValueError: If substitution matrix is not provided.
    :raises ImportError: If Biopython is not available.
    """
    if not BIOPYTHON_AVAILABLE:
        raise ImportError("Biopython is required for using the aligner. Please install Biopython to use this functionality!")

    if mode not in ["global", "local"]:
        raise ValueError("mode must be one of 'global' or 'local'!")

    aligner = PairwiseAligner()
    aligner.mode = mode
    aligner.open_internal_insertion_score = open_internal_insertion_score
    aligner.extend_internal_insertion_score = extend_internal_insertion_score
    aligner.open_left_insertion_score = open_left_insertion_score
    aligner.extend_left_insertion_score = extend_left_insertion_score
    aligner.open_right_insertion_score = open_right_insertion_score
    aligner.extend_right_insertion_score = extend_right_insertion_score
    aligner.open_internal_deletion_score = open_internal_deletion_score
    aligner.extend_internal_deletion_score = extend_internal_deletion_score
    aligner.open_left_deletion_score = open_left_deletion_score
    aligner.extend_left_deletion_score = extend_left_deletion_score
    aligner.open_right_deletion_score = open_right_deletion_score
    aligner.extend_right_deletion_score = extend_right_deletion_score
    aligner.wildcard = None
    aligner.substitution_matrix = substitution_matrix

    return aligner
