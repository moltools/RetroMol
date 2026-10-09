"""Module defining the Result data class."""

from dataclasses import dataclass
from typing import Any

from retromol.model.submission import Submission
from retromol.model.reaction_graph import ReactionGraph
from retromol.model.readout import LinearReadout
from retromol.chem.mol import smiles_to_mol
from retromol.chem.tagging import get_tags_mol


@dataclass(frozen=True)
class Result:
    """
    Represents a RetroMol parsing result.

    :var submission: Submission: The original submission associated with this result.
    :var reaction_graph: ReactionGraph: The reaction graph generated from retrosynthetic analysis.
    :var linear_readout: LinearReadout: The linear readout representation of the reaction graph.
    """

    submission: Submission
    reaction_graph: ReactionGraph
    linear_readout: LinearReadout

    def __str__(self) -> str:
        """
        String representation of the Result.
        
        :return: String representation of the Result.
        """
        return f"Result(submission={self.submission}, reaction_graph={self.reaction_graph}, linear_readout={self.linear_readout})"
    
    def calculate_coverage(self) -> float:
        """
        Calculate the fraction of input heavy atoms in identified assembly monomers.

        :return: Coverage score as a float.
        """
        input_heavy_atoms = smiles_to_mol(self.submission.smiles).GetNumHeavyAtoms()
        if not input_heavy_atoms:
            return 0.0

        root_heavy_tags = {
            atom.GetIsotope()
            for atom in self.submission.mol.GetAtoms()
            if atom.GetAtomicNum() > 1 and atom.GetIsotope() != 0
        }
        identified_tags: set[int] = set()
        for node in self.linear_readout.assembly_graph.monomer_nodes():
            if node.is_identified:
                identified_tags.update(get_tags_mol(node.mol))

        return len(identified_tags.intersection(root_heavy_tags)) / input_heavy_atoms

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the Result to a dictionary.

        :return: Dictionary representation of the Result.
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

        :param data: Dictionary representation of the Result.
        :return: Result object.
        """
        submission = Submission.from_dict(data["submission"])
        reaction_graph = ReactionGraph.from_dict(data["reaction_graph"])
        linear_readout = LinearReadout.from_dict(data["linear_readout"])

        return cls(
            submission=submission,
            reaction_graph=reaction_graph,
            linear_readout=linear_readout,
        )
