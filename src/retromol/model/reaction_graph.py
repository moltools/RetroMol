"""Data structures for representing reaction application graphs."""

import logging
import json
from collections import defaultdict, deque
from dataclasses import dataclass, field, replace
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

    :var enc: str: Unique encoding of the molecule.
    :var mol: Mol: The molecule object.
    :var smiles: str: SMILES representation of the molecule.
    :var identity: MolIdentity | None: Identification information if identified.
    :var identified: bool | None: Whether the node has been checked for identification.
    """

    enc: str
    mol: Mol
    smiles: str
    identity: MolIdentity | None = None
    identified: bool | None = None  # None=unknown, False=checked-no, True=checked-yes
    pruned: bool = False  # uninformative descendants were collapsed here

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
        
        :return: str: String representation of the MolNode.
        """
        id_name = self.identity.name if self.identity else None
        return f"MolNode(enc={self.enc}, id={id_name})"

    def identify(self, rules: list[MatchingRule], match_stereochemistry: bool = False) -> MolIdentity | None:
        """
        Identify the molecule node using the provided matching rules.

        :param rules: The matching rules to apply.
        :param match_stereochemistry: Whether to consider stereochemistry in matching.
        :return: The identity if matched, else None.
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

        :return: Dictionary representation of the MolNode.
        """
        return {
            "enc": self.enc,
            "tagged_smiles": mol_to_smiles(self.mol, include_tags=True),
            "smiles": self.smiles,
            "identity": self.identity.to_dict() if self.identity else None,
            "identified": self.identified,
            "pruned": self.pruned,
        }
    
    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "MolNode":
        """
        Deserialize a MolNode from a dictionary.

        :param data: Dictionary representation of the MolNode.
        :return: MolNode object.
        """
        identity = MolIdentity.from_dict(data["identity"]) if data["identity"] else None

        node = cls(
            enc=data["enc"],
            mol=smiles_to_mol(data["tagged_smiles"]),
            smiles=data["smiles"],
            identity=identity,
            identified=data["identified"],
            pruned=data.get("pruned", False),
        )
        return node


@dataclass(frozen=True)
class ReactionStep:
    """
    Edge payload: desribes one application event:
    - uncontested: multiple rules applied as one step
    - contested: exactly one rule applied

    :var kind: StepKind: 'uncontested' or 'contested'
    :var names: Tuple[str, ...]: reaction rule IDs (human-facing).
    :var rule_ids: Tuple[str, ...]: optional numeric IDs (stable internal).
    """

    kind: StepKind
    names: tuple[str, ...]  # reaction rule IDs (human-facing)
    rule_ids: tuple[str, ...] = ()  # optional numeric IDs (stable internal)

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the ReactionStep to a dictionary.

        :return: Dictionary representation of the ReactionStep.
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

        :param data: Dictionary representation of the ReactionStep.
        :return: ReactionStep object.
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

    :var src: str: Encoding of source molecule node.
    :var dsts: Tuple[str, ...]: encodings of child molecule nodes.
    :var step: ReactionStep: details of the reaction application.
    """

    src: str
    dsts: tuple[str, ...]
    step: ReactionStep

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the RxnEdge to a dictionary.

        :return: Dictionary representation of the RxnEdge.
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

        :param data: Dictionary representation of the RxnEdge.
        :return: RxnEdge object.
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
        
        :return: Mapping of encodings to identified MolNodes.
        """
        return {enc: node for enc, node in self.nodes.items() if node.is_identified}

    def __str__(self) -> str:
        """
        Return a string representation of the ReactionGraph.
        
        :return: String representation.
        """
        return f"ReactionGraph(num_nodes={len(self.nodes)}, num_edges={len(self.edges)})"

    def add_node(self, mol: Mol) -> str:
        """
        Add a molecule node to the graph if not already present.
        
        :param mol: Molecule to add.
        :param keep_stereo_smiles: Whether to keep stereochemistry in SMILES.
        :return: Encoding of the molecule node.
        """
        enc = encode_mol(mol)
        if enc not in self.nodes:
            self.nodes[enc] = MolNode(enc=enc, mol=Mol(mol), smiles=mol_to_smiles(mol, include_tags=False))
            self.out_edges.setdefault(enc, [])

        return enc
    
    def add_edge(self, src_enc: str, child_mols: Iterable[Mol], step: ReactionStep) -> tuple[str, ...]:
        """
        Add a reaction edge to the graph.

        :param src_enc: Encoding of the source molecule node.
        :param child_mols: Iterable of child molecule nodes.
        :param step: ReactionStep describing the reaction.
        :return: Tuple of encodings of the child molecule nodes.
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

        :param identified_only: Whether to include only identified nodes.
        :return: List of MolNode objects that are leaves.
        """
        leaves: list[MolNode] = []

        for enc, node in self.nodes.items():
            # No outgoing edges -> leaf
            if not self.out_edges.get(enc):
                if identified_only and not node.is_identified:
                    continue
                leaves.append(node)

        return leaves

    def prune_unidentified(self, root_enc: str) -> "ReactionGraph":
        """Collapse wholly unidentified subtrees to unresolved boundary nodes.

        Keep every explored alternative leading to an identification, together
        with all its sibling products. Identified intermediates remain selectable.
        Direct self-dependent reactions cannot form finite decompositions and
        are discarded. Longer cycles are handled when selecting assemblies.
        """
        if root_enc not in self.nodes:
            raise ValueError(f"Root {root_enc!r} is absent from the reaction graph")
        incoming: dict[str, set[str]] = defaultdict(set)
        outgoing: dict[str, list[RxnEdge]] = defaultdict(list)
        for edge in self.edges:
            if edge.src in edge.dsts:
                continue
            outgoing[edge.src].append(edge)
            for child in edge.dsts:
                incoming[child].add(edge.src)

        productive = set(self.identified_nodes)
        queue = deque(productive)
        while queue:
            for parent in incoming[queue.popleft()]:
                if parent not in productive:
                    productive.add(parent)
                    queue.append(parent)

        graph = ReactionGraph()
        queue = deque([root_enc])
        while queue:
            enc = queue.popleft()
            if enc in graph.nodes:
                continue
            node = self.nodes[enc]
            # A wholly unproductive reaction alternative adds no information.
            # Retain unresolved products only when they are siblings required
            # by an informative reaction, rather than creating cheap dead-end
            # alternatives by truncating those reactions to unknown leaves.
            edges = [edge for edge in outgoing[enc] if any(d in productive for d in edge.dsts)]
            collapsed = bool(self.out_edges.get(enc)) and not edges
            graph.nodes[enc] = replace(node, pruned=True) if collapsed else node
            graph.out_edges[enc] = []
            for edge in edges:
                graph.out_edges[enc].append(len(graph.edges))
                graph.edges.append(edge)
                queue.extend(edge.dsts)
        return graph

    def to_dict(self) -> dict[str, Any]:
        """Compact graph: indexed nodes/edges and interned identities/steps.

        Atom-tagged SMILES preserve structure, stereo, and original-atom mapping.
        Encodings, untagged SMILES and adjacency are reconstructed on load.
        Node indices follow sorted encodings, also used by Result.root.
        """
        indices = {enc: i for i, enc in enumerate(sorted(self.nodes))}
        identities: list[dict[str, Any]] = []
        identity_indices: dict[str, int] = {}
        nodes = []
        for enc in indices:
            node = self.nodes[enc]
            record: dict[str, Any] = {"smiles": mol_to_smiles(node.mol, include_tags=True)}
            if node.is_identified:
                if node.identity is None:
                    raise ValueError(f"Identified node {enc!r} has no identity")
                identity = node.identity.matched_rule.to_dict()
                key = json.dumps(identity, sort_keys=True, separators=(",", ":"))
                if key not in identity_indices:
                    identity_indices[key] = len(identities)
                    identities.append(identity)
                record["identity"] = identity_indices[key]
            elif not node.is_checked:
                record["checked"] = False
            if node.pruned:
                record["pruned"] = True
            nodes.append(record)

        steps: list[dict[str, Any]] = []
        step_indices: dict[ReactionStep, int] = {}
        edges = []
        for edge in self.edges:
            if edge.step not in step_indices:
                step_indices[edge.step] = len(steps)
                steps.append(edge.step.to_dict())
            edges.append({"src": indices[edge.src], "dsts": [indices[d] for d in edge.dsts],
                          "step": step_indices[edge.step]})
        return {"nodes": nodes, "identities": identities, "steps": steps, "edges": edges}

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ReactionGraph":
        """Load the compact graph format; rebuild encodings and adjacency."""
        identities = [MolIdentity(MatchingRule.from_dict(d)) for d in data["identities"]]
        steps = [ReactionStep.from_dict(d) for d in data["steps"]]
        graph = cls()
        encodings = []
        for record in data["nodes"]:
            mol = smiles_to_mol(record["smiles"])
            enc = encode_mol(mol)
            if enc in graph.nodes:
                raise ValueError("Duplicate molecule in reaction graph nodes")
            identity = cls._at(identities, record["identity"]) if "identity" in record else None
            graph.nodes[enc] = MolNode(
                enc, mol, mol_to_smiles(mol), identity,
                True if identity else (False if record.get("checked", True) else None),
                pruned=record.get("pruned", False),
            )
            graph.out_edges[enc] = []
            encodings.append(enc)
        for record in data["edges"]:
            edge = RxnEdge(cls._at(encodings, record["src"]),
                           tuple(cls._at(encodings, d) for d in record["dsts"]),
                           cls._at(steps, record["step"]))
            graph.out_edges[edge.src].append(len(graph.edges))
            graph.edges.append(edge)
        return graph

    @staticmethod
    def _at(items: list, index: int):
        if type(index) is not int or not 0 <= index < len(items):
            raise ValueError(f"Invalid graph reference: {index!r}")
        return items[index]
