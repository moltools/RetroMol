"""Module defining the Result data class."""

from dataclasses import dataclass, asdict
from typing import Any


@dataclass(frozen=True)
class Result:
    """
    Represents a RetroMol parsing result.
    """

    ...

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the Result to a dictionary.

        :return: dictionary representation of the Result
        """
        return asdict(self)
