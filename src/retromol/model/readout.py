"""Data structures for representing readouts from RetroMol parsing results."""

from dataclasses import dataclass
from typing import Literal

from retromol.model.reaction_graph import MolNode, ReactionGraph
from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.rules import MatchingRule
from retromol.chem.mol import encode_mol
from retromol.chem.tagging import get_tags_mol


ReadoutMode = Literal["leaf_identified", "first_identified"]


@dataclass(frozen=True)
class LinearReadout:
    """
    A linear readout representation of a RetroMol parsing result.
    """

    assembly_graph: AssemblyGraph
    paths: list[list[MolNode]]

    def __str__(self) -> str:
        """
        Return a string representation of the LinearReadout.

        :return: str: string representation
        """
        return f"LinearReadout(assembly_graph_nodes={self.assembly_graph.g.number_of_nodes()}, assembly_graph_edges={self.assembly_graph.g.number_of_edges()}, num_paths={len(self.paths)})"

    @classmethod
    def from_reaction_graph(
        cls,
        root_enc: str,
        reaction_graph: ReactionGraph,
        exclude_identities: list[MatchingRule] | None = None,
        include_identities: list[MatchingRule] | None = None,
    ) -> "LinearReadout":
        """
        Create a LinearReadout from a Result object.

        :param root_enc: encoding of the root molecule
        :param reaction_graph: ReactionGraph object
        :param exclude_identities: list of matching rules to exclude identities (not used here)
        :param include_identities: list of matching rules to include identities (not used here)
        :return: LinearReadout instance
        :raises ValueError: if root_enc not found in reaction graph nodes
        """
        exclude_identities = exclude_identities or []
        include_identities = include_identities  # keep None meaning "no whitelist"

        # Convert identities to their IDs for easier checking
        exclude_identities = set([r.id for r in exclude_identities])
        if include_identities is not None:
            include_identities = set([r.id for r in include_identities])

        g = reaction_graph
        if root_enc not in g.nodes:
            raise ValueError(f"root_enc {root_enc} not found in reaction graph nodes")
        
        # Use root_enc to get root mol
        root = g.nodes[root_enc].mol

        # Create assembly graph of monomers; first collect nodes to include
        collected = g.get_leaf_nodes(identified_only=False)
        a = AssemblyGraph.build(root_mol=root, monomers=collected, include_unassigned=True)

        # Break bonds between monomers that are not backbone-related bonds (i.e., keep C-C and C-N bonds only)
        f = a.filtered_by_root_bond_elements(allow_pairs={frozenset(("C", "C")), frozenset(("C", "N"))}, drop_isolated=False)

        # Get individual connected components from the assembly graph and extract longest paths (allow to visit each edge only once)
        hs = f.connected_components()

        paths: list[list[MolNode]] = []
        for h in hs:
            path = h.longest_path()
            paths.append(path)

        return cls(assembly_graph=a, paths=paths)
    
    def to_dict(self) -> dict:
        """
        Serialize the LinearReadout to a dictionary.

        :return: dict representation of LinearReadout
        """
        return {
            "assembly_graph": self.assembly_graph.to_dict(),
            "paths": [[node.to_dict() for node in path] for path in self.paths],
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> "LinearReadout":
        """
        Deserialize a LinearReadout from a dictionary.

        :param data: dict representation of LinearReadout
        :return: LinearReadout instance
        """
        assembly_graph = AssemblyGraph.from_dict(data["assembly_graph"])
        paths = [[MolNode.from_dict(node_data) for node_data in path_data] for path_data in data["paths"]]
        return cls(assembly_graph=assembly_graph, paths=paths)
