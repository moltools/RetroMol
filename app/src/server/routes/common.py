# -*- coding: utf-8 -*-

"""This module contains the common settings for the application."""

import os
import typing as ty
from enum import Enum, auto


NEO4J_URI = os.environ.get("NEO4J_URI", "bolt://localhost:7687")
NEO4J_USER = os.environ.get("NEO4J_USER", "neo4j")
NEO4J_PASSWORD = os.environ.get("NEO4J_PASSWORD", "password")


class ResponseStatus(Enum):
    """The status of a response.

    :cvar Success: The response was successful.
    :cvar Warning: The response was successful, but with a warning.
    :cvar Failure: The response failed.
    """
    Success = auto()
    Warning = auto()
    Failure = auto()

    def __str__(self) -> str:
        return self.name.lower()


class ResponseData:
    """A class for creating a response object.
    
    :param status: The status of the response.
    :type status: ResponseStatus
    :param payload: The payload of the response.
    :type payload: dict
    :param message: The message of the response.
    """
    def __init__(
        self,
        status: ty.Optional[ResponseStatus] = None,
        payload: ty.Optional[dict] = None,
        message: ty.Optional[str] = None
    ) -> None:
        # Cut off message if it is too long.
        if message is not None and len(message) > 80:
            message = message[:80] + "..."

        self.status = status if status is not None else ResponseStatus.Warning
        self.payload = payload if payload is not None else dict()
        self.message = message if message is not None else "No message provided!"

    def jsonify(self) -> ty.Dict[str, ty.Any]:
        """Return the response as a dictionary.
        
        :return: The response as a dictionary.
        :rtype: ty.Dict[str, ty.Any]
        """
        return dict(
            status=str(self.status),
            payload=self.payload,
            message=self.message
        )
    
def fail(message: ty.Optional[str] = None) -> ty.Dict[str, ty.Any]: 
    """Return a failed response.
    
    :param message: The message of the response.
    :type message: str
    :return: The response.
    :rtype: ty.Dict[str, ty.Any]
    """
    return ResponseData(status=ResponseStatus.Failure, message=message).jsonify()

def warning(message: ty.Optional[str] = None) -> ty.Dict[str, ty.Any]:
    """Return a warning response.
    
    :param message: The message of the response.
    :type message: str
    :return: The response.
    :rtype: ty.Dict[str, ty.Any]
    """
    return ResponseData(status=ResponseStatus.Warning, message=message).jsonify()

def success(
    message: ty.Optional[str] = None,
    payload: ty.Optional[ty.Dict[str, ty.Any]] = None
) -> ty.Dict[str, ty.Any]:
    """Return a successful response.
    
    :param message: The message of the response.
    :type message: str
    :param payload: The payload of the response.
    :type payload: ty.Dict[str, ty.Any]
    :return: The response.
    :rtype: ty.Dict[str, ty.Any]
    """
    return ResponseData(status=ResponseStatus.Success, payload=payload, message=message).jsonify()
