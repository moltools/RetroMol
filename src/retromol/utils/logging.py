"""Utility functions for logging."""

from __future__ import annotations

import logging
import sys


PACKAGE_LOGGER = "retromol"

STANDARD_FMT = "%(asctime)s | %(levelname)s | %(name)s | %(message)s"
STANDARD_DATEFMT = "%Y-%m-%d %H:%M:%S"


def setup_logging(
    level: str | int = "INFO",
    *,
    fmt: str = STANDARD_FMT,
    datefmt: str = STANDARD_DATEFMT,
    stream: None | int | str | object = None,
) -> None:
    """
    Set up logging for the retromol package.

    :param level: Log level for console output.
    :param fmt: Log message format.
    :param datefmt: Date format for log messages.
    :param stream: Output stream for console logs; defaults to sys.stderr.
    .. note:: Safe to call multiple times; library code should not call this function;
        it is intended for use by applications using the library.
    """
    if stream is None:
        stream = sys.stderr

    if isinstance(level, str):
        level = level.upper()

    root = logging.getLogger()
    root.setLevel(level)

    handler = logging.StreamHandler(stream)
    handler.setLevel(level)
    handler.setFormatter(logging.Formatter(fmt=fmt, datefmt=datefmt))

    # Avoid duplicate handlers if called repeatedly (common in notebooks)
    # Keep it simple: remove existing handlers created by previous setup calls.
    root.handlers = [handler]

    # Make sure package logger propagates to root
    logging.getLogger(PACKAGE_LOGGER).propagate = True


def add_file_handler(
    logfile: str,
    *,
    level: str | int = "DEBUG",
    fmt: str = STANDARD_FMT,
    datefmt: str = STANDARD_DATEFMT,
) -> None:
    """
    Add a file handler to the root logger.

    :param logfile: Path to log file.
    :param level: Log level for file output.
    .. note:: Intended to be called after setup_logging(); safe to call multiple times
        for the same logfile.
    """
    if isinstance(level, str):
        level = level.upper()

    root = logging.getLogger()

    # Prevent duplicate file handlers for the same path
    for h in root.handlers:
        if isinstance(h, logging.FileHandler) and h.baseFilename == logfile:
            return

    fh = logging.FileHandler(logfile)
    fh.setLevel(level)
    fh.setFormatter(logging.Formatter(fmt=fmt, datefmt=datefmt))

    root.addHandler(fh)
