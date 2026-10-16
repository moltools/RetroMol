"""Data structures for representing readouts from RetroMol parsing results."""

from dataclasses import dataclass
from typing import Callable, Literal, TypeVar

from retromol.model.reaction_graph import MolNode, ReactionGraph
from retromol.model.assembly_graph import AssemblyGraph
from retromol_fingerprint.fingerprint import TOKEN_LINK


ReadoutMode = Literal["leaf_identified", "first_identified"]

T = TypeVar("T")


def merge_named_paths(paths: list[list[T]], key: Callable[[T], str], link_item: T) -> list[T]:
    """
    Merge candidate primary sequence paths into a single sequence, joined by
    `link_item`, so a compound/BGC with more than one disconnected path (e.g. a
    branched/cyclic assembly, or a separately represented tailoring moiety)
    still gets exactly one primary
    sequence to store/query/align, instead of one entry per path.

    Longer paths sort first; paths of equal length are ordered lexicographically
    by their own token sequence (via `key`), so the merge is deterministic and
    reproducible regardless of the order `paths` was produced in.

    :param paths: candidate paths, each a list of items (e.g. names, or (name, tags) pairs)
    :param key: extracts the string used for sorting/tie-breaking from one item
    :param link_item: the item inserted between two merged paths (e.g. TOKEN_LINK,
        or a (TOKEN_LINK, []) pair matching the shape of the other items)
    :return: one flat list: every path's items in sorted order, `link_item` between them
    """
    ordered = sorted(paths, key=lambda items: (-len(items), [key(item) for item in items]))

    merged: list[T] = []
    for i, items in enumerate(ordered):
        if i > 0:
            merged.append(link_item)
        merged.extend(items)

    return merged


def split_named_paths(sequence: list[T], key: Callable[[T], str], link_token: str = TOKEN_LINK) -> list[list[T]]:
    """
    Inverse of `merge_named_paths`: split a merged sequence back apart into its
    constituent chains, wherever an item's `key` equals `link_token`. The link
    items themselves are dropped, never included in a returned chain -- so
    `merge_named_paths(split_named_paths(seq, key, tok), key, link_item)` round-trips
    `seq` (given a `link_item` whose `key` is `tok`).

    :param sequence: a merged sequence, as returned by `merge_named_paths`
    :param key: extracts the string used to test for the link token from one item
    :param link_token: the token marking a join point (default: TOKEN_LINK)
    :return: the chains, in their original relative order, empty chains dropped
    """
    chains: list[list[T]] = [[]]
    for item in sequence:
        if key(item) == link_token:
            chains.append([])
        else:
            chains[-1].append(item)

    return [chain for chain in chains if chain]


@dataclass(frozen=True)
class LinearReadout:
    """
    A linear readout representation of a RetroMol parsing result.

    :cvar assembly_graph: AssemblyGraph representing the assembly of monomers.
    :cvar paths: Disjoint paths of identified MolNodes. Unresolved regions break
        paths and remain only in the source assembly graph.
    """

    assembly_graph: AssemblyGraph
    paths: list[list[MolNode]]

    def __post_init__(self) -> None:
        """Keep direct construction and older serialized paths identified-only."""
        paths: list[list[MolNode]] = []
        for path in self.paths:
            part: list[MolNode] = []
            for node in path:
                if node.is_identified and node.identity is not None:
                    part.append(node)
                elif part:
                    paths.append(part)
                    part = []
            if part:
                paths.append(part)
        object.__setattr__(self, "paths", paths)

    def __str__(self) -> str:
        """
        Return a string representation of the LinearReadout.

        :return: str: String representation.
        """
        return f"LinearReadout(assembly_graph_nodes={self.assembly_graph.g.number_of_nodes()}, assembly_graph_edges={self.assembly_graph.g.number_of_edges()}, num_paths={len(self.paths)})"

    @classmethod
    def from_reaction_graph(
        cls,
        root_enc: str,
        reaction_graph: ReactionGraph,
        identified_only: bool = False,
        include_unassigned: bool = True,
    ) -> "LinearReadout":
        """
        Create a LinearReadout from one selected reaction route.

        An explored graph can contain mutually exclusive, overlapping leaves.
        Select a route first, or use from_assembly_graph with a sampled assembly.

        :param root_enc: Encoding of the root molecule.
        :param reaction_graph: ReactionGraph object.
        :param identified_only: Whether to only include identified monomers in
            the source assembly graph. Readout paths are always identified-only.
        :param include_unassigned: Whether to include an unassigned node in the assembly graph.
        :return: LinearReadout instance.
        :raises ValueError: If root_enc not found in reaction graph nodes.
        """
        g = reaction_graph
        if root_enc not in g.nodes:
            raise ValueError(f"The root_enc {root_enc} not found in reaction graph nodes!")
        
        # Use root_enc to get root mol
        root = g.nodes[root_enc].mol

        # Create assembly graph of monomers; first collect nodes to include
        collected = g.get_leaf_nodes(identified_only=identified_only)
        a = AssemblyGraph.build(root_mol=root, monomers=collected, include_unassigned=include_unassigned)
        return cls.from_assembly_graph(a)

    @classmethod
    def from_assembly_graph(cls, assembly_graph: AssemblyGraph) -> "LinearReadout":
        """Partition a sampled assembly into deterministic, disjoint paths.

        Matching-rule ``sequence_role`` separates tailoring units from other
        units. Connections between units with the same explicit role survive
        regardless of bond elements (e.g. a connected sugar chain). Otherwise,
        the C-C/C-N heuristic applies to non-tailoring units.

        Repeatedly extract a long simple path and process its remaining branches
        so every identified MolNode occurs exactly once. Unidentified and
        unassigned nodes are excluded before traversal, breaking paths at gaps.
        The input graph is unchanged; use its coverage to measure completeness.
        Rings are opened deterministically, without assigning a biological start
        or direction; cyclic path selection is heuristic.
        """
        g = assembly_graph.g.copy()
        if g.has_node(assembly_graph.unassigned):
            g.remove_node(assembly_graph.unassigned)

        roles = {}
        for node_id, data in list(g.nodes(data=True)):
            node = data.get("molnode")
            if node is None:
                raise ValueError(f"AssemblyGraph node {node_id!r} has no 'molnode' attached!")
            if not node.is_identified or node.identity is None:
                g.remove_node(node_id)
                continue
            roles[node_id] = node.identity.matched_rule.sequence_role

        backbone_pairs = {frozenset(("C", "C")), frozenset(("C", "N"))}
        for u, v, data in list(g.edges(data=True)):
            if roles[u] is not None and roles[u] == roles[v]:
                continue
            if "tailoring" in (roles[u], roles[v]) or not any(
                frozenset((bond.a1_symbol, bond.a2_symbol)) in backbone_pairs
                for bond in data["bonds"]
            ):
                g.remove_edge(u, v)

        paths: list[list[MolNode]] = []
        pending = AssemblyGraph(g, unassigned=assembly_graph.unassigned).connected_components()
        while pending:
            h = pending.pop()
            path = h.longest_path()
            paths.append(path)
            # Graph node IDs need not be identical to MolNode.enc.
            used = {id(node) for node in path}
            h.g.remove_nodes_from([n for n, data in h.g.nodes(data=True) if id(data["molnode"]) in used])
            pending.extend(h.connected_components())

        paths.sort(key=lambda path: (-len(path), [_node_name(node) for node in path], [node.enc for node in path]))
        return cls(assembly_graph=assembly_graph, paths=paths)

    def primary_sequence_parts(self) -> list[list[str]]:
        """Named backbone/tailoring paths, longest first with lexical tie-breaks.

        Each identified assembly monomer occurs once, including repeated
        identities. Unresolved and unassigned regions are omitted; assembly
        coverage measures completeness. Sublists do not encode attachment
        points, omitted crosslinks, or ring closures; retain ``assembly_graph``
        when those are needed. Untagged paths are not biological role labels.
        """
        return sorted(
            [[_node_name(node) for node in path] for path in self.paths],
            key=lambda path: (-len(path), path),
        )

    def primary_sequence(self) -> list[str]:
        """
        Compatibility flat view of ``primary_sequence_parts()``, joined by
        TOKEN_LINK. Use the nested view to keep separate backbone, branch and
        tailoring paths. Only identified monomers are included.

        :return: the single merged primary sequence
        """
        return merge_named_paths(self.primary_sequence_parts(), key=lambda n: n, link_item=TOKEN_LINK)

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

        :param data: Dictionary representation of LinearReadout.
        :return: LinearReadout instance.
        """
        assembly_graph = AssemblyGraph.from_dict(data["assembly_graph"])
        paths = [[MolNode.from_dict(node_data) for node_data in path_data] for path_data in data["paths"]]
        return cls(assembly_graph=assembly_graph, paths=paths)


def _node_name(node: MolNode) -> str:
    return node.identity.name


def assembly_to_primary_sequence(assembly_graph: AssemblyGraph) -> list[list[str]]:
    """Return linear parts of one assembly, including backbone and tailoring.

    Uses matching-rule sequence roles with a C-C/C-N fallback and includes every
    identified monomer once. Unresolved and unassigned regions are excluded;
    use ``result.calculate_coverage(assembly_graph)`` for completeness. Parts are ordered
    longest first, then lexically. Ordering is structural, not biosynthetic;
    branch attachments and ring closures remain in the input assembly graph.

    For the corresponding MolNodes (including atom tags and identities), use
    ``LinearReadout.from_assembly_graph(assembly_graph).paths``.
    """
    return LinearReadout.from_assembly_graph(assembly_graph).primary_sequence_parts()
