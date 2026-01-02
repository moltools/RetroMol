"""Data structures for representing reaction application graphs."""

import logging
from dataclasses import dataclass, field
from typing import Iterable, Literal

from retromol.chem.mol import Mol, encode_mol
from retromol.model.identity import MolIdentity
from retromol.model.rules import MatchingRule
from retromol.chem.matching import match_mol

log = logging.getLogger(__name__)


StepKind = Literal["uncontested", "contested"]


@dataclass(frozen=True)
class MolNode:
    """
    A molecule node in the processing graph.

    :var enc: str: unique encoding of the molecule (canonical SMILES with isotopic tags)
    :var mol: Mol: the molecule object
    :var identity: MolIdentity | None: identification information if identified
    :var identified: bool | None: whether the node has been checked for identification
    """

    enc: str
    mol: Mol
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


@dataclass(frozen=True)
class ReactionStep:
    """
    Edge payload: desribes one application event.
    - uncontested: multiple rules applied as one step
    - contested: exactly one rule applied

    :var kind: StepKind: 'uncontested' or 'contested'
    :var names: tuple[str, ...]: reaction rule IDs (human-facing)
    :var rule_ids: tuple[int, ...]: optional numeric IDs (stable internal)
    """

    kind: StepKind
    names: tuple[str, ...]  # reaction rule IDs (human-facing)
    rule_ids: tuple[int, ...] = ()  # optional numeric IDs (stable internal)


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


@dataclass
class ReactionGraph:
    """
    Simple directed hypergraph:
    - nodes: enc -> MolNode
    - edges: list of RxnEdge
    - out_edges: adjacency index for fast traversal
    """

    nodes: dict[int, MolNode] = field(default_factory=dict)
    edges: list[RxnEdge] = field(default_factory=list)
    out_edges: dict[int, list[int]] = field(default_factory=dict)  # enc -> indices into edges

    @property
    def identified_nodes(self) -> dict[int, MolNode]:
        """
        Return only identified nodes.
        
        :return: dict[int, MolNode]: mapping of encodings to identified MolNodes
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
        :return: encoding of the molecule node
        """
        enc = encode_mol(mol)
        if enc not in self.nodes:
            self.nodes[enc] = MolNode(enc=enc, mol=Mol(mol))
            self.out_edges.setdefault(enc, [])

        return enc
    
    def add_edge(self, src_enc: int, child_mols: Iterable[Mol], step: ReactionStep) -> tuple[int, ...]:
        """
        Add a reaction edge to the graph.

        :param src_enc: encoding of the source molecule node
        :param child_mols: iterable of child molecule nodes
        :param step: ReactionStep describing the reaction
        :return: tuple of encodings of the child molecule nodes
        """
        dst_encs: list[int] = []
        for m in child_mols:
            dst_encs.append(self.add_node(m))

        edge = RxnEdge(src=src_enc, dsts=tuple(dst_encs), step=step)
        self.edges.append(edge)
        self.out_edges.setdefault(src_enc, []).append(len(self.edges) - 1)

        return tuple(dst_encs)
    