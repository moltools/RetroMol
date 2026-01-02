"""Utilities for function timeouts using SIGALRM."""

import signal
from collections.abc import Callable, Generator
from typing import Any, Generic, ParamSpec, TypeVar


P = ParamSpec("P")
T = TypeVar("T")


class FunctionTimeoutError(Exception):
    """
    Custom exception for function timeout.
    """

    pass


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
