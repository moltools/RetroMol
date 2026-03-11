"""Module contains funcionalities for comparing fingerprints."""

import numpy as np


def tanimoto_binary(fp1: np.ndarray, fp2: np.ndarray) -> float:
    """
    Compute the Tanimoto similarity between two binary fingerprints.

    :param fp1: First fingerprint (binary).
    :param fp2: Second fingerprint (binary).
    :return: Tanimoto similarity score between 0 and 1.
    """
    a = np.bitwise_and(fp1, fp2).sum()
    b = np.bitwise_or(fp1, fp2).sum()
    return float(a / b) if b > 0 else 0.0


def cosine_counted(fp1: np.ndarray, fp2: np.ndarray) -> float:
    """
    Compute the cosine similarity between two counted fingerprints.

    :param fp1: First fingerprint (counted).
    :param fp2: Second fingerprint (counted).
    :return: Cosine similarity score between 0 and 1.
    """
    dot = np.dot(fp1, fp2)
    norm = np.linalg.norm(fp1) * np.linalg.norm(fp2)
    return float(dot / norm) if norm > 0 else 0.0
