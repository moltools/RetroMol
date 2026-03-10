"""Fingerprint similarity metrics."""

import numpy as np
from numpy.typing import NDArray


def calculate_cosine_similarity(fp1: NDArray[np.int8], fp2: NDArray[np.int8]) -> float:
    """
    Cosine similarity for fingerprints.

    :param fp1: first fingerprint (1D array)
    :param fp2: second fingerprint (1D array)
    :return: cosine similarity in [0, 1]
    """
    a = np.asarray(fp1)
    b = np.asarray(fp2)

    # Ensure 1D
    a = a.ravel()
    b = b.ravel()
    if a.shape != b.shape:
        raise ValueError(f"Different lengths: {a.shape} vs {b.shape}")

    # Upcast to float to avoid integer overflow and match sklearn
    a = a.astype(np.float64, copy=False)
    b = b.astype(np.float64, copy=False)

    # Compute cosine
    dot = float(np.dot(a, b))
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))

    if na == 0.0 or nb == 0.0:
        return 0.0

    return dot / (na * nb)


def calculate_tanimoto_similarity(fp1: NDArray[np.int8], fp2: NDArray[np.int8]) -> float:
    """
    Tanimoto similarity for molecular fingerprints (binary or count-based).

    :param fp1: first fingerprint (1D array)
    :param fp2: second fingerprint (1D array)
    :return: Tanimoto similarity in [0, 1]
    """
    a = np.asarray(fp1)
    b = np.asarray(fp2)

    # Ensure 1D
    a = a.ravel()
    b = b.ravel()
    if a.shape != b.shape:
        raise ValueError(f"Different lengths: {a.shape} vs {b.shape}")

    # Upcast to float to prevent overflow and ensure precision
    a = a.astype(np.float64, copy=False)
    b = b.astype(np.float64, copy=False)

    # Dot product = intersection term
    ab = float(np.dot(a, b))
    aa = float(np.dot(a, a))
    bb = float(np.dot(b, b))

    denom = aa + bb - ab
    if denom == 0.0:
        return 0.0

    return ab / denom
