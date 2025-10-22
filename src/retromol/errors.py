"""This module defines custom exceptions for the RetroMol package."""


class FunctionTimeoutError(Exception):
    """Custom exception for function timeout."""

    pass


class MotifGraphNodeWithoutAttributesError(Exception):
    """Custom error for when a motif graph node is created without attributes."""

    def __init__(self, message: str) -> None:
        """
        Initialize the error with a message.

        :param message: the error message
        """
        super().__init__(message)
        self.message = message

    def __str__(self) -> str:
        """
        Return the string representation of the error.

        :return: the error message
        """
        return f"MotifGraphNodeWithoutAttributesError: {self.message}"
