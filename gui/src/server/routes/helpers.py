"""Module providing helper functions for endpoints."""

import uuid
from typing import Any

import numpy as np

BITS_PER_HEX_DIGIT = 4


def get_unique_identifier() -> str:
    """
    Generate a unique identifier string.

    :return: unique identifier as a string
    """
    return str(uuid.uuid4())


def bits_to_hex(bits: np.ndarray, n_bits: int = 512) -> str:
    """
    Convert a numpy array of bits (0/1 ints) into a hexadecimal string representation.

    :param bits: numpy array of shape (n_bits,) or (1, n_bits) with values 0 or 1
    :param n_bits: expected number of bits (default is 512)
    :return: hexadecimal string representation
    :raises ValueError: if input array shape is incorrect or contains invalid values
    """
    arr = np.asarray(bits, dtype=np.int8).reshape(-1)

    if arr.size != n_bits:
        raise ValueError(f"Input array must have shape ({n_bits},) or (1, {n_bits})")

    if n_bits % BITS_PER_HEX_DIGIT != 0:
        raise ValueError("Number of bits must be a multiple of 4 in order to convert to hexadecimal")

    # Guard that values are actually 0/1
    if not np.isin(arr, (0, 1)).all():
        raise ValueError("Input array must only contain 0 and 1 values")

    bitstring = "".join("1" if b else "0" for b in arr)  # n_bits characters
    hex_len = n_bits // BITS_PER_HEX_DIGIT
    hexstr = format(int(bitstring, 2), f"0{hex_len}x")

    return hexstr


def hex_to_bits(hexstr: str, n_bits: int = 512) -> np.ndarray:
    """
    Convert a hexadecimal string representation into a numpy array of bits (0/1 ints).

    :param hexstr: hexadecimal string representation
    :param n_bits: expected number of bits (default is 512)
    :return: numpy array of shape (n_bits,) with values 0 or 1
    :raises ValueError: if input string is invalid or does not match expected bit length
    """
    hexstr = hexstr.strip().lower()

    # Basic hex validation
    if not hexstr:
        raise ValueError("Input hexadecimal string is empty")
    if any(c not in "0123456789abcdef" for c in hexstr):
        raise ValueError("Input string contains non-hexadecimal characters")

    inferred_bits = len(hexstr) * BITS_PER_HEX_DIGIT
    if inferred_bits != n_bits:
        raise ValueError(
            f"Input hexadecimal string must represent {n_bits} bits (length {n_bits // BITS_PER_HEX_DIGIT})")

    bit_int = int(hexstr, 16)
    bitstring = bin(bit_int)[2:].zfill(n_bits)  # binary string of length n_bits

    return np.fromiter((1 if b == "1" else 0 for b in bitstring), dtype=np.int8, count=n_bits)


def kmerize_sequence(sequence: list[Any], k: int) -> list[list[Any]]:
    """
    Generate k-mers from a given sequence (forward and backward).

    :param sequence: list of elements (e.g., amino acids)
    :param k: length of each k-mer
    :return: list of k-mer strings
    """
    kmers = []
    seq_length = len(sequence)

    # Forward k-mers
    for i in range(seq_length - k + 1):
        kmer = sequence[i:i + k]
        kmers.append(kmer)

    # Backward k-mers
    for i in range(seq_length - k, -1, -1):
        kmer = sequence[i:i + k]
        kmers.append(kmer)

    return kmers