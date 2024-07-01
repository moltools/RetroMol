# -*- coding: utf-8 -*-

"""This module contains helper functions for RetroMol."""

import errno
import functools
import os
import signal
import typing as ty

TIMEOUT_ERROR_MSG = os.strerror(errno.ETIME)


class RetroMolTimeoutError(Exception):
    """Custom exception to raise when a function exceeds the specified time."""

    pass


def timeout(seconds: int = 5, error_message: str = TIMEOUT_ERROR_MSG) -> ty.Callable:
    """Raise a RetroMolTimeoutError when runtime exceeds the specified time.

    :param seconds: The number of seconds before the function times out.
    :type seconds: int
    :param error_message: The error message to display when the function times
        out.
    :type error_message: str
    :return: Decorator.
    :rtype: ty.Callable
    """

    def decorator(func: ty.Callable) -> ty.Callable:
        """Decorate the function."""

        def _handle_timeout(signum, frame):
            raise RetroMolTimeoutError(error_message)

        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            """Wrap the function."""
            signal.signal(signal.SIGALRM, _handle_timeout)
            signal.alarm(seconds)
            try:
                result = func(*args, **kwargs)
            finally:
                signal.alarm(0)
            return result

        return wrapper

    return decorator
