# -*- coding: utf-8 -*-

"""Helper functions for the NPKG package."""

import logging
import os
import typing as ty


def validate_path(path: str) -> None:
    """Check if a path exists.

    :param path: The path to check.
    :type path: str
    :raises FileNotFoundError: If the path does not exist.
    """
    logger = logging.getLogger(__name__)

    if os.path.exists(path):
        msg = f"Path {path} exists."
        logger.info(msg)

    else:
        msg = f"Path {path} does not exist."
        logger.error(msg)
        raise FileNotFoundError(msg)


def safe_get(
    dictionary: dict, keys: ty.List[str], default: ty.Optional[ty.Any] = None
) -> ty.Optional[ty.Any]:
    """Safely get a value from a dictionary.

    :param dictionary: The dictionary to search.
    :type dictionary: dict
    :param keys: The keys to search for.
    :type keys: ty.List[str]
    :param default: The default value if the key is not found.
    :type default: ty.Optional[ty.Any]
    :return: The value if found, otherwise None.
    :rtype: ty.Any
    """
    for key in keys:
        try:
            dictionary = dictionary[key]
        except (KeyError, TypeError):
            return default

    return dictionary
