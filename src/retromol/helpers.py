"""This module provides helper functions for the RetroMol package."""

import hashlib
import json
import signal
from collections.abc import Callable, Generator
from typing import Any, Generic, ParamSpec, TypeVar

import ijson

from retromol.errors import FunctionTimeoutError

P = ParamSpec("P")
T = TypeVar("T")


def sha256_hex(s: str) -> str:
    """
    Compute the SHA-256 hash of the input string `s` and return its hexadecimal representation.
    If `s` is None, treat it as an empty string.

    :param s: input string to hash
    :return: hexadecimal representation of the SHA-256 hash
    """
    return hashlib.sha256((s or "").encode("utf-8")).hexdigest()


def blake64_hex(s: str) -> str:
    """
    Compute the BLAKE2b hash of the input string `s` and return the first 16 characters

    :param s: input string to hash
    :return: 16-character hex string (64 bits)
    """
    return hashlib.blake2b((s or "").encode("utf-8"), digest_size=8).hexdigest()


def _timeout_handler(signum: int, frame: Any) -> None:
    """
    Signal handler that raises our custom timeout exception.

    :param signum: the signal number
    :param frame: the current stack frame
    """
    raise FunctionTimeoutError("function timed out")


class _TimeoutWrapper(Generic[P, T]):
    """
    A callable wrapper that runs `func` under a SIGALRM timeout.
    Because this class is defined at module scope, instances are picklable.
    """

    __slots__ = ("func", "seconds")

    def __init__(self, func: Callable[P, T], seconds: int) -> None:
        """
        Initialize the _TimeoutWrapper.

        :param func: the function to wrap
        :param seconds: the timeout duration in seconds
        """
        self.func = func
        self.seconds = seconds

    def __call__(self, *args: P.args, **kwargs: P.kwargs) -> T:
        """
        Call the wrapped function with a timeout.

        :param args: positional arguments for the wrapped function
        :param kwargs: keyword arguments for the wrapped function
        :return: the result of the wrapped function
        :raises FunctionTimeoutError: if the function call times out
        """
        original_handler = signal.getsignal(signal.SIGALRM)
        signal.signal(signal.SIGALRM, _timeout_handler)
        signal.alarm(self.seconds)

        try:
            return self.func(*args, **kwargs)
        finally:
            # Cancel any pending alarm and restore the original handler
            signal.alarm(0)
            signal.signal(signal.SIGALRM, original_handler)


def timeout_decorator(seconds: int = 30) -> Callable[[Callable[P, T]], _TimeoutWrapper[P, T]]:
    """
    Returns a decorator that wraps a function in a _TimeoutWrapper instance.

    :param seconds: The timeout duration in seconds
    :return: A decorator that applies the timeout to a function
    .. note:: because _TimeoutWrapper is at module level, it is picklable
    """

    def decorate(func: Callable[P, T]) -> _TimeoutWrapper[P, T]:
        return _TimeoutWrapper(func, seconds)

    return decorate


def iter_json(path: str, jsonl: bool = False) -> Generator[Any, None, None]:
    """
    Stream items from a JSON array or a JSON Lines (JSONL) file.

    :param path: path to the JSON or JSONL file
    :param jsonl: if True, treat the file as JSONL (one JSON object per line). If False, assume a single JSON array
    :yield: parsed JSON objects
    """
    with open(path, "rb") as f:
        if jsonl:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                yield json.loads(line)
        else:
            yield from ijson.items(f, "item")
