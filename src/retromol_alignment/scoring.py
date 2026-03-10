"""Scoring functionalities for sequence alignment."""

import numpy as np
import pandas as pd

from retromol_alignment.aligner import substitution_matrices


def create_scoring_matrix(df: pd.DataFrame) -> substitution_matrices.Array:
    """
    Create a scoring matrix from a pandas DataFrame.

    :param df: pandas DataFrame containing the scoring matrix, where the index and columns are the same and represent the alphabet, and the values are the scoring values
    :return: scoring matrix as a substitution_matrices.Array object
    """
    if df.shape[0] != df.shape[1]:
        raise ValueError("Substitution matrix must be square (same number of rows and columns).")

    alphabet = tuple(df.columns)
    data = df.to_numpy(dtype=np.float64)

    sm = substitution_matrices.Array(alphabet, 2, data, np.float64)

    return sm
