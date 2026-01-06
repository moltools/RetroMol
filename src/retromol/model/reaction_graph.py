"""Data structures for representing reaction application graphs."""

import logging
from dataclasses import dataclass, field
from typing import Any, Iterable, Literal

from retromol.chem.mol import Mol, encode_mol, mol_to_smiles, smiles_to_mol
from retromol.model.identity import MolIdentity
from retromol.model.rules import MatchingRule
from retromol.chem.matching import match_mol

log = logging.getLogger(__name__)


StepKind = Literal["uncontested", "contested"]


@dataclass(frozen=True)
class MolNode:
    """
    A molecule node in the processing graph.

    :var enc: str: unique encoding of the molecule
    :var mol: Mol: the molecule object
    :var smiles: str: SMILES representation of the molecule
    :var identity: MolIdentity | None: identification information if identified
    :var identified: bool | None: whether the node has been checked for identification
    """

    enc: str
    mol: Mol
    smiles: str
    identity: MolIdentity | None = None
    identified: bool | None = None  # None=unknown, False=checked-no, True=checked-yes

    @property
    def is_checked(self) -> bool:
        return self.identified is not None
    
    @property
    def is_identified(self) -> bool:
        return self.identified is True
    
    @property
    def is_unidentified_checked(self) -> bool:
        return self.identified is False
    
    def __str__(self) -> str:
        """
        Return a string representation of the MolNode.
        
        :return: str: string representation
        """
        id_name = self.identity.name if self.identity else None
        return f"MolNode(enc={self.enc}, id={id_name})"

    def identify(self, rules: list[MatchingRule], match_stereochemistry: bool = False) -> MolIdentity | None:
        """
        Identify the molecule node using the provided matching rules.

        :param rules: list[MatchingRule]: the matching rules to apply
        :param match_stereochemistry: bool: whether to consider stereochemistry in matching
        :return: MolIdentity | None: the identity if matched, else None
        """
        if self.is_checked:
            return self.identity  # identity is present only if identified=True

        if identity := match_mol(self.mol, rules, match_stereochemistry):
            object.__setattr__(self, "identity", identity)
            object.__setattr__(self, "identified", True)
            return identity
        
        object.__setattr__(self, "identity", None)
        object.__setattr__(self, "identified", False)
        return None
    
    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the MolNode to a dictionary.

        :return: dictionary representation of the MolNode
        """
        return {
            "enc": self.enc,
            "tagged_smiles": mol_to_smiles(self.mol, include_tags=True),
            "smiles": self.smiles,
            "identity": self.identity.to_dict() if self.identity else None,
            "identified": self.identified,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "MolNode":
        """
        Deserialize a MolNode from a dictionary.

        :param data: dictionary representation of the MolNode
        :return: MolNode object
        """
        identity = MolIdentity.from_dict(data["identity"]) if data["identity"] else None

        node = cls(
            enc=data["enc"],
            mol=smiles_to_mol(data["tagged_smiles"]),
            smiles=data["smiles"],
            identity=identity,
            identified=data["identified"],
        )
        return node


@dataclass(frozen=True)
class ReactionStep:
    """
    Edge payload: desribes one application event.
    - uncontested: multiple rules applied as one step
    - contested: exactly one rule applied

    :var kind: StepKind: 'uncontested' or 'contested'
    :var names: tuple[str, ...]: reaction rule IDs (human-facing)
    :var rule_ids: tuple[str, ...]: optional numeric IDs (stable internal)
    """

    kind: StepKind
    names: tuple[str, ...]  # reaction rule IDs (human-facing)
    rule_ids: tuple[str, ...] = ()  # optional numeric IDs (stable internal)

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the ReactionStep to a dictionary.

        :return: dictionary representation of the ReactionStep
        """
        return {
            "kind": self.kind,
            "names": self.names,
            "rule_ids": self.rule_ids,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ReactionStep":
        """
        Deserialize a ReactionStep from a dictionary.

        :param data: dictionary representation of the ReactionStep
        :return: ReactionStep object
        """
        step = cls(
            kind=data["kind"],
            names=tuple(data["names"]),
            rule_ids=tuple(data.get("rule_ids", ())),
        )
        return step


@dataclass
class RxnEdge:
    """
    Directed hyper-edge parent -> children, labeled by ReactionStep.

    :var src: int: encoding of source molecule node
    :var dsts: tuple[int, ...]: encodings of child molecule nodes
    :var step: ReactionStep: details of the reaction application
    """

    src: int
    dsts: tuple[int, ...]
    step: ReactionStep

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the RxnEdge to a dictionary.

        :return: dictionary representation of the RxnEdge
        """
        return {
            "src": self.src,
            "dsts": self.dsts,
            "step": self.step.to_dict(),
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "RxnEdge":
        """
        Deserialize a RxnEdge from a dictionary.

        :param data: dictionary representation of the RxnEdge
        :return: RxnEdge object
        """
        edge = cls(
            src=data["src"],
            dsts=tuple(data["dsts"]),
            step=ReactionStep.from_dict(data["step"]),
        )
        return edge


@dataclass
class ReactionGraph:
    """
    Simple directed hypergraph:
    - nodes: enc -> MolNode
    - edges: list of RxnEdge
    - out_edges: adjacency index for fast traversal
    """

    nodes: dict[str, MolNode] = field(default_factory=dict)
    edges: list[RxnEdge] = field(default_factory=list)
    out_edges: dict[str, list[int]] = field(default_factory=dict)  # enc -> indices into edges

    @property
    def identified_nodes(self) -> dict[str, MolNode]:
        """
        Return only identified nodes.
        
        :return: dict[str, MolNode]: mapping of encodings to identified MolNodes
        """
        return {enc: node for enc, node in self.nodes.items() if node.is_identified}

    def __str__(self) -> str:
        """
        Return a string representation of the ReactionGraph.
        
        :return: str: string representation
        """
        return f"ReactionGraph(num_nodes={len(self.nodes)}, num_edges={len(self.edges)})"

    def add_node(self, mol: Mol) -> int:
        """
        Add a molecule node to the graph if not already present.
        
        :param mol: molecule to add
        :param keep_stereo_smiles: whether to keep stereochemistry in SMILES
        :return: encoding of the molecule node
        """
        enc = encode_mol(mol)
        if enc not in self.nodes:
            self.nodes[enc] = MolNode(enc=enc, mol=Mol(mol), smiles=mol_to_smiles(mol, include_tags=False))
            self.out_edges.setdefault(enc, [])

        return enc
    
    def add_edge(self, src_enc: str, child_mols: Iterable[Mol], step: ReactionStep) -> tuple[str, ...]:
        """
        Add a reaction edge to the graph.

        :param src_enc: encoding of the source molecule node
        :param child_mols: iterable of child molecule nodes
        :param step: ReactionStep describing the reaction
        :return: tuple of encodings of the child molecule nodes
        """
        dst_encs: list[str] = []
        for m in child_mols:
            dst_encs.append(self.add_node(m))

        edge = RxnEdge(src=src_enc, dsts=tuple(dst_encs), step=step)
        self.edges.append(edge)
        self.out_edges.setdefault(src_enc, []).append(len(self.edges) - 1)

        return tuple(dst_encs)
    
    def get_leaf_nodes(self, identified_only: bool = True) -> list[MolNode]:
        """
        Get all leaf nodes (nodes with no outgoing edges).

        :param identified_only: whether to include only identified nodes
        :return: list of MolNode objects that are leaves
        """
        leaves: list[MolNode] = []

        for enc, node in self.nodes.items():
            # No outgoing edges -> leaf
            if not self.out_edges.get(enc):
                if identified_only and not node.is_identified:
                    continue
                leaves.append(node)

        return leaves

    
    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the ReactionGraph to a dictionary.

        :return: dictionary representation of the ReactionGraph
        """
        return {
            "nodes": {enc: node.to_dict() for enc, node in self.nodes.items()},
            "edges": [edge.to_dict() for edge in self.edges],
            "out_edges": {enc: indices for enc, indices in self.out_edges.items()},
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ReactionGraph":
        """
        Deserialize a ReactionGraph from a dictionary.

        :param data: dictionary representation of the ReactionGraph
        :return: ReactionGraph object
        """
        reaction_graph = cls(
            nodes={enc: MolNode.from_dict(node_data) for enc, node_data in data["nodes"].items()},
            edges=[RxnEdge.from_dict(edge_data) for edge_data in data["edges"]],
            out_edges={enc: indices for enc, indices in data["out_edges"].items()},
        )
        return reaction_graph
