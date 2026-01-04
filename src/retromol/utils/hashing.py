"""Module for hashing utilities in RetroMol."""

import hashlib


def sha256_hex(s: str) -> str:
    """
    Compute the SHA-256 hash of string and return its hexadecimal representation.
    
    :param s: input string to hash
    :return: hexadecimal representation of the SHA-256 hash
    .. note:: None is treated as an empty string
    """
    return hashlib.sha256((s or "").encode("utf-8")).hexdigest()


def blake64_hex(s: str) -> str:
    """
    Compute the BLAKE2b hash of string and return the first 16 characters.
    
    :param s: input string to hash
    :return: first 16 characters of the BLAKE2b hash in hexadecimal representation
    .. note:: None is treated as an empty string
    """
    return hashlib.blake2b((s or "").encode("utf-8"), digest_size=8).hexdigest()
