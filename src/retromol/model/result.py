"""Module defining the Result data class."""

from dataclasses import dataclass
from typing import Any

from retromol.model.submission import Submission
from retromol.model.reaction_graph import ReactionGraph
from retromol.model.readout import LinearReadout
from retromol.chem.tagging import get_tags_mol


@dataclass(frozen=True)
class Result:
    """
    Represents a RetroMol parsing result.

    :var submission: Submission: the original submission associated with this result
    :var reaction_graph: ReactionGraph: the reaction graph generated from retrosynthetic analysis
    :var linear_readout: LinearReadout: the linear readout representation of the reaction graph
    """

    submission: Submission
    reaction_graph: ReactionGraph
    linear_readout: LinearReadout

    def __str__(self) -> str:
        """
        String representation of the Result.
        
        :return: string representation of the Result
        """
        return f"Result(submission={self.submission}, reaction_graph={self.reaction_graph}, linear_readout={self.linear_readout})"
    
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
            "linear_readout": self.linear_readout.to_dict(),
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
        linear_readout = LinearReadout.from_dict(data["linear_readout"])

        return cls(
            submission=submission,
            reaction_graph=reaction_graph,
            linear_readout=linear_readout,
        )
