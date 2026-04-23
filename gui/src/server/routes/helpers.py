"""Module providing helper functions for endpoints."""

import uuid


def get_unique_identifier() -> str:
    """
    Generate a unique identifier string.

    :return: unique identifier as a string
    """
    return str(uuid.uuid4())
