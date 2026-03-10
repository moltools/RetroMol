"""Multiple Sequence Alignment (MSA) module."""

import logging

import numpy as np

from retromol_alignment.aligner import PairwiseAligner
from retromol_alignment.pairwise import T, Converter, _pairwise_alignment


log = logging.getLogger(__name__)


def calculate_msa(
    aligner: PairwiseAligner,
    to_align: list[list[T]],
    converter: Converter,
    center_star: int | None = None,
    can_reverse: list[bool] | None = None
) -> tuple[list[tuple[float, bool, list[T | None]]], list[int]]:
    """
    Calculate a multiple sequence alignment (MSA) for the given sequences.

    :param aligner: PairwiseAligner object to use for pairwise alignments
    :param to_align: list of sequences to align, where each sequence is a list of items
    :param converter: Converter object to convert items in sequences to integers for alignment and back
    :param center_star: index of the sequence to use as the center for the star alignment (default: None, which means to
        use the first sequence)
    :param can_reverse: list of booleans indicating whether each sequence can be reversed for alignment (default: None,
        which means none of the sequences can be reversed)
    :return: tuple containing a list of tuples, where each tuple contains the alignment score, a boolean indicating if
        the sequence was reversed,and the aligned sequence (with None for gaps),and a list of indices representing the
        order of sequences in the MSA
    :raises ValueError: if the length of can_reverse does not match the number of sequences to align
    """
    if not to_align:
        return []
    
    if can_reverse is None:
        can_reverse = [False] * len(to_align)

    if len(can_reverse) != len(to_align):
        raise ValueError(f"Length of can_reverse ({len(can_reverse)}) must match the number of sequences to align ({len(to_align)}).")

    # Failsafe: set aligner mode to global
    alignment_mode = aligner.mode
    if alignment_mode != "global":
        log.warning(f"Aligner mode is '{alignment_mode}', but 'global' is required for MSA. Setting aligner mode to 'global'.")
        aligner.mode = "global"

    # Convert sequences into int arrays based on substitution matrix alphabet
    int_seqs: list[np.ndarray[np.int32]] = [converter.to_int_array(seq) for seq in to_align]

    # Create pairwise similarity matrix; use reverse orientation if it scores better
    sims = np.zeros((len(int_seqs), len(int_seqs)), dtype=np.float32)
    best_reversed = np.zeros((len(int_seqs), len(int_seqs)), dtype=bool)

    for i, int_seq1 in enumerate(int_seqs):
        for j, int_seq2 in enumerate(int_seqs[i+1:], start=i+1):
            score = aligner.score(int_seq1, int_seq2)
            is_reversed = False

            if can_reverse[j]:
                rev_score = aligner.score(int_seq1, int_seq2[::-1].copy())
                if rev_score > score:
                    score = rev_score
                    is_reversed = True
            
            # We only calculate the lower triangle of the similarity matrix
            if i >= j:
                continue

            sims[i, j] = score
            sims[j, i] = score
            best_reversed[i, j] = is_reversed
            best_reversed[j, i] = is_reversed

    # Mask similarity matrix to ignore self-similarity (diagonal)
    masked_sims = sims.copy()
    np.fill_diagonal(masked_sims, -np.inf)

    # Find center star sequence if not given, ignore the diagonal (self-similarity)
    if center_star is None:
        center_ind = int(np.argmax(masked_sims.sum(axis=0)))
    else:
        center_ind = int(center_star)

    # Sort all by descending similarity to center star
    sorted_indices = np.argsort(masked_sims[center_ind])[::-1]
    other_inds = [i for i in sorted_indices if i != center_ind]

    # Align every sequence to the center star sequence
    scores = []
    reversed_flags = []
    msa = np.array([int_seqs[center_ind]], dtype=np.int32)

    for other_ind in other_inds:
        # Align other sequence to center star sequence
        t = msa[0]
        q = int_seqs[other_ind]

        # Use the same orientation choice as already determined
        if best_reversed[center_ind, other_ind]:
            q = q[::-1].copy()
            is_reversed = True
        else:
            is_reversed = False

        s, t_a, q_a = _pairwise_alignment(aligner, t, q, converter.insert_repr)

        # Insert gaps in remainder sequences, if any
        remainder = msa[1:]
        if remainder.shape[0] > 0:
            insert_inds = [i for i, x in enumerate(t_a) if x == converter.insert_repr]
            for insert_ind in insert_inds:
                remainder = np.insert(remainder, insert_ind, converter.insert_repr, axis=1)
            msa = np.vstack([t_a, remainder, q_a])
        else:
            msa = np.array([t_a, q_a], dtype=np.int32)
        
        scores.append(s)
        reversed_flags.append(is_reversed)

        # Rename insert_repr to gap_repr
        msa[msa == converter.insert_repr] = converter.gap_repr

    # Double check if there are any gap-only columns, delete those
    gap_columns = np.all(msa == converter.gap_repr, axis=0)
    msa = msa[:, ~gap_columns]

    # Convert back to original sequences with None for gaps
    aligned_seqs: list[tuple[float, list[T | None]]] = []
    for i, int_seq in enumerate(msa):
        aligned_seq = converter.from_int_array(int_seq)

        if i == 0:
            # Center star sequence, score is self-alignment score
            score = sims[center_ind, center_ind]
        else:
            # Score is the pairwise alignment score to the center star sequence
            score = scores[i - 1]

        is_rev = reversed_flags[i - 1] if i > 0 else False

        aligned_seqs.append((score, is_rev, aligned_seq))

    new_order = [center_ind] + other_inds

    return aligned_seqs, new_order
