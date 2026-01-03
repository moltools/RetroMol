"""Module defining the Result data class."""

from dataclasses import dataclass, asdict
from typing import Any

from retromol.model.submission import Submission
from retromol.model.reaction_graph import ReactionGraph
from retromol.chem.tagging import get_tags_mol


@dataclass(frozen=True)
class Result:
    """
    Represents a RetroMol parsing result.
    """

    submission: Submission
    reaction_graph: ReactionGraph

    def __str__(self) -> str:
        """
        String representation of the Result.
        
        :return: string representation of the Result
        """
        return f"Result(submission={self.submission}, reaction_graph={self.reaction_graph})"
    
    def calculate_coverage(self) -> float:
        """
        Calculate coverage score for result.
        
        :return: coverage score as a float
        """
        # Collect all unique tags from identified nodes
        identified_tags = set()
        for node in self.reaction_graph.identified_nodes.values():
            identified_tags.update(get_tags_mol(node.mol))

        # Get all unique tags from the root
        root_tags = set(get_tags_mol(self.submission.mol))

        # Calculate coverage: proportion of root tags identified
        if root_tags:
            coverage = len(identified_tags.intersection(root_tags)) / len(root_tags)
            return coverage

        return 0.0

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the Result to a dictionary.

        :return: dictionary representation of the Result
        """
        return {
            "submission": self.submission.to_dict(),
            "reaction_graph": self.reaction_graph.to_dict(),
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Result":
        """
        Deserialize a Result from a dictionary.

        :param data: dictionary representation of the Result
        :return: Result object
        """
        submission = Submission.from_dict(data["submission"])
        reaction_graph = ReactionGraph.from_dict(data["reaction_graph"])

        return cls(
            submission=submission,
            reaction_graph=reaction_graph,
        )
