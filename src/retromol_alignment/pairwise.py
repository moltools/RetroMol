"""Pairwise sequence alignment module."""

from dataclasses import dataclass
from typing import TypeVar, Callable

import numpy as np
from numpy.typing import NDArray

from retromol_alignment.aligner import PairwiseAligner


T = TypeVar("T")



@dataclass(frozen=True)
class Converter:
    """
    Converter for sequences to be aligned.

    :cvar to_identifier: function to convert items in sequences to integers for alignment
    :cvar from_identifier: function to convert integers back to items in sequences after alignment
    :cvar gap_repr: integer representation of a gap in the alignment (default: -1)
    :cvar insert_repr: integer representation of an insertion in the alignment (default: -
    """

    to_identifier: Callable[[T], np.int32]
    from_identifier: Callable[[np.int32], T]

    gap_repr: np.int32 = np.int32(-1)
    insert_repr: np.int32 = np.int32(-2)

    def to_int_array(self, sequence: list[T]) -> NDArray[np.int32]:
        """
        Convert a sequence to a numpy array of integers for alignment.

        :param sequence: list of items in the sequence
        :return: numpy array of integers representing the sequence
        """
        return np.array([self.to_identifier(item) for item in sequence], dtype=np.int32)
    
    def from_int_array(self, int_array: NDArray[np.int32]) -> list[T | None]:
        """
        Convert a numpy array of integers back to a sequence of items after alignment.

        :param int_array: numpy array of integers representing the aligned sequence (with gap_repr for gaps)
        :return: list of items in the aligned sequence (with None for gaps)
        """
        return [self.from_identifier(item) if item != self.gap_repr else None for item in int_array]


def _pairwise_alignment(
    aligner: PairwiseAligner,
    t: NDArray[np.int32],
    q: NDArray[np.int32],
    gap_repr: np.int32 = np.int32(-1),
) -> tuple[float, NDArray[np.int32], NDArray[np.int32]]:
    """
    Align two sequences and return the alignment result.

    :param aligner: PairwiseAligner object to use for alignment
    :param t: first sequence as a numpy array of integers
    :param q: second sequence as a numpy array of integers
    :param gap_repr: integer representation of a gap in the alignment (default: -1)
    :return: tuple containing the alignment score, aligned first sequence, and aligned second sequence
    :raises ValueError: if the alignment coordinates are invalid
    """
    alignments = aligner.align(seqA=t, seqB=q)

    # Pick first alignment
    alignment = alignments[0]
    score = alignments[0].score

    t_a: list[np.int32] = []
    q_a: list[np.int32] = []
    for i in range(alignment.coordinates.shape[1] - 1):
        a = alignment.coordinates[0][i:i+2]
        b = alignment.coordinates[1][i:i+2]
        len_a = a[1] - a[0]
        len_b = b[1] - b[0]
        if len_a == len_b:
            t_a.extend(t[a[0]:a[1]])
            q_a.extend(q[b[0]:b[1]])
        elif len_a == 0:
            t_a.extend([gap_repr] * len_b)
            q_a.extend(q[b[0]:b[1]])
        elif len_b == 0:
            t_a.extend(t[a[0]:a[1]])
            q_a.extend([gap_repr] * len_a)
        else:
            raise ValueError("Invalid alignment coordinates")
        
    return score, np.array(t_a, dtype=np.int32), np.array(q_a, dtype=np.int32)


def align(
    aligner: PairwiseAligner,
    t: list[T],
    q: list[T],
    converter: Converter,
) -> tuple[float, list[T | None], list[T | None]]:
    """
    Align two sequences and return the alignment result.

    :param aligner: PairwiseAligner object to use for alignment
    :param t: first sequence as a list of items
    :param q: second sequence as a list of items
    :param converter: Converter object to convert items in sequences to integers for alignment and back
    :return: tuple containing the alignment score, aligned first sequence, and aligned second sequence
    """
    t_int = np.array([converter.to_identifier(item) for item in t], dtype=np.int32)
    q_int = np.array([converter.to_identifier(item) for item in q], dtype=np.int32)

    s, t_a, q_a = _pairwise_alignment(aligner=aligner, t=t_int, q=q_int, gap_repr=converter.gap_repr)

    t_a_converted = [converter.from_identifier(item) if item != converter.gap_repr else None for item in t_a]
    q_a_converted = [converter.from_identifier(item) if item != converter.gap_repr else None for item in q_a]

    return s, t_a_converted, q_a_converted
