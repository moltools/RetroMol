"""Data structures for representing readouts from RetroMol parsing results."""

from dataclasses import dataclass
from typing import Callable, Literal, TypeVar

import networkx as nx

from retromol.chem.mol import smarts_to_mol
from retromol.model.reaction_graph import MolNode, ReactionGraph
from retromol.model.assembly_graph import AssemblyGraph
from retromol_fingerprint.fingerprint import TOKEN_LINK


ReadoutMode = Literal["leaf_identified", "first_identified"]

T = TypeVar("T")

_CARBOXYL = smarts_to_mol("[CX3](=[OX1])-[OX2H1,OX1-]")
_UNIT_ENTRY_CAP = smarts_to_mol("[#6]-[0S]-[0OH]")


def _carboxyl_exit_tags(node: MolNode, include_added_caps: bool = True) -> set[int]:
    """Original atoms at a fragment's carboxyl-side continuation points.

    A retained carboxyl carbon is the donor in a cleaved peptide/ester link.
    When the entire carboxyl is an added cap (tag zero), its original carbon
    neighbor is the continuation point instead. This uses the same fragment
    representation for peptide, PK, and hybrid paths, without rematching units.
    """
    exits: set[int] = set()
    for carbon_idx, oxo_idx, hydroxyl_idx in node.mol.GetSubstructMatches(_CARBOXYL):
        carbon = node.mol.GetAtomWithIdx(carbon_idx)
        if carbon.GetIsotope():
            exits.add(carbon.GetIsotope())
        elif include_added_caps and all(node.mol.GetAtomWithIdx(i).GetIsotope() == 0
                                        for i in (oxo_idx, hydroxyl_idx)):
            exits.update(atom.GetIsotope() for atom in carbon.GetNeighbors()
                         if atom.GetAtomicNum() == 6 and atom.GetIsotope())
    return exits


def _remove_side_links(g: nx.Graph) -> None:
    """Omit side attachments that compete with connected backbone entries.

    The artificial S-OH cap marks the original atom receiving the preceding
    unit. If that entry is connected, an additional acyl link to the same unit's
    O/S is a side attachment, whether it closes a ring or joins separate chains.
    Being a graph bridge does not make that attachment a backbone continuation.
    An extra ether link to the already connected incoming carbon likewise does
    not replace the connection from the preceding unit's carboxyl side.
    Esters without this competing entry evidence retain their connectivity.
    This mutates only the readout's graph copy, never the source assembly.
    """
    nodes = {key: data["molnode"] for key, data in g.nodes(data=True)
             if data["molnode"].identity.matched_rule.sequence_role != "tailoring"}
    entries = {
        key: {node.mol.GetAtomWithIdx(match[0]).GetIsotope()
              for match in node.mol.GetSubstructMatches(_UNIT_ENTRY_CAP)
              if node.mol.GetAtomWithIdx(match[0]).GetIsotope()}
        for key, node in nodes.items()
    }
    if not any(entries.values()):
        return
    exits = {key: _carboxyl_exit_tags(node) for key, node in nodes.items()}
    retained = {key: _carboxyl_exit_tags(node, include_added_caps=False) for key, node in nodes.items()}
    connected_entries = {key: set() for key in nodes}
    for u, v, data in g.edges(data=True):
        if u not in nodes or v not in nodes:
            continue
        for bond in data["bonds"]:
            for left, right in ((bond.a1_tag, bond.a2_tag), (bond.a2_tag, bond.a1_tag)):
                if left in exits[u] and right in entries[v]:
                    connected_entries[v].add(right)
                if left in entries[u] and right in exits[v]:
                    connected_entries[u].add(left)
    if not any(connected_entries.values()):
        return

    for u, v, data in list(g.edges(data=True)):
        if u not in nodes or v not in nodes:
            continue
        kept = []
        for bond in data["bonds"]:
            side_attachment = False
            for carbon, hetero, symbol in (
                (bond.a1_tag, bond.a2_tag, bond.a2_symbol),
                (bond.a2_tag, bond.a1_tag, bond.a1_symbol),
            ):
                if bond.bond_type == "SINGLE" and symbol in ("O", "S") and (
                    carbon in retained[u] and hetero in g.nodes[v]["tags"] and connected_entries[v]
                    or carbon in retained[v] and hetero in g.nodes[u]["tags"] and connected_entries[u]
                    or carbon in connected_entries[u] and hetero in g.nodes[v]["tags"]
                    or carbon in connected_entries[v] and hetero in g.nodes[u]["tags"]
                ):
                    side_attachment = True
            if not side_attachment:
                kept.append(bond)
        if len(kept) == len(data["bonds"]):
            continue
        if kept:
            g[u][v].update(bonds=kept, n_bonds=len(kept))
        else:
            g.remove_edge(u, v)


def _readout_graph(assembly: AssemblyGraph) -> nx.Graph:
    """Shared connectivity for both path selection and path orientation."""
    g = assembly.g.copy()
    if g.has_node(assembly.unassigned):
        g.remove_node(assembly.unassigned)
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
    _remove_side_links(g)
    return g


def _orient_paths(assembly: AssemblyGraph, paths: list[list[MolNode]]) -> list[list[MolNode]]:
    """Orient paths toward the carboxyl end when connection evidence agrees.

    Root bond tags identify the actual connection; merely finding a free acid
    (including an Asp/Glu side-chain acid) is not a direction signal. Missing or
    contradictory connection evidence can still be resolved by a unique retained
    terminal carboxyl. Otherwise the existing deterministic orientation is kept.
    Closed rings keep their selected start, a convention rather than a terminus.
    """
    g = _readout_graph(assembly)
    by_enc = {data["molnode"].enc: node_id for node_id, data in g.nodes(data=True)}
    exits = {node.enc: _carboxyl_exit_tags(node) for path in paths for node in path}
    oriented = []
    for path in paths:
        if len(path) < 2 or any(node.identity.matched_rule.sequence_role == "tailoring" for node in path):
            oriented.append(path)
            continue
        ids = [by_enc.get(node.enc) for node in path]
        if any(node_id is None for node_id in ids):
            oriented.append(path)
            continue

        closed = len(path) > 2 and g.has_edge(ids[-1], ids[0])
        steps = list(zip(range(len(path) - 1), range(1, len(path))))
        if closed:
            steps.append((len(path) - 1, 0))
        directions: set[int] = set()
        for i, j in steps:
            edge = g.get_edge_data(ids[i], ids[j])
            if edge is None:
                continue
            left_tags = g.nodes[ids[i]]["tags"]
            right_tags = g.nodes[ids[j]]["tags"]
            for bond in edge["bonds"]:
                # Root bond begin/end order is unrelated to sequence direction.
                for left, right in ((bond.a1_tag, bond.a2_tag), (bond.a2_tag, bond.a1_tag)):
                    if left in left_tags and right in right_tags:
                        if left in exits[path[i].enc]:
                            directions.add(1)
                        if right in exits[path[j].enc]:
                            directions.add(-1)
        reverse = directions == {-1}
        if len(directions) != 1 and not closed:
            terminal_ends = []
            for i, j in ((0, 1), (len(path) - 1, len(path) - 2)):
                retained = _carboxyl_exit_tags(path[i], include_added_caps=False)
                edge = g.get_edge_data(ids[i], ids[j], {})
                connected = {tag for bond in edge.get("bonds", []) for tag in (bond.a1_tag, bond.a2_tag)}
                # A side-chain acid is not a terminal anchor when the same unit
                # donates through another retained carboxyl to its neighbor.
                if retained and not retained.intersection(connected):
                    terminal_ends.append(i)
            if len(terminal_ends) == 1:
                reverse = terminal_ends[0] == 0
        if reverse:
            path = list(reversed(path))
            if closed:
                # Preserve the chosen ring start while reversing its traversal.
                path = path[-1:] + path[:-1]
        oriented.append(path)
    return oriented


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
        """Normalize identified paths and their direction, including saved paths."""
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
        paths = _orient_paths(self.assembly_graph, paths)
        paths.sort(key=lambda path: (-len(path), [_node_name(node) for node in path], [node.enc for node in path]))
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
        Ester closures competing with an already connected cap-defined unit
        entry are omitted from traversal before selecting paths.
        Carboxyl-side connections orient consistent paths N/starter-to-C end.
        Ambiguous paths retain deterministic direction. Rings have a conventional
        start, not a biological terminus; cyclic path selection is heuristic.
        """
        g = _readout_graph(assembly_graph)
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

        return cls(assembly_graph=assembly_graph, paths=paths)

    def primary_sequence_parts(self) -> list[list[str]]:
        """Named backbone/tailoring paths, longest first with lexical tie-breaks.

        Each identified assembly monomer occurs once, including repeated
        identities. Unresolved and unassigned regions are omitted; assembly
        coverage measures completeness. Sublists do not encode attachment
        points, omitted crosslinks, or ring closures; retain ``assembly_graph``
        when those are needed. Untagged paths are not biological role labels.
        Consistent carboxyl-side connection evidence puts the N/starter end first.
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
    longest first, then lexically. Within each path, consistent carboxyl-side
    connections orient N/starter-to-C end; ambiguous paths retain deterministic
    order. This is not proof of biosynthesis. Branch attachments and ring closures
    remain in the input assembly graph.

    For the corresponding MolNodes (including atom tags and identities), use
    ``LinearReadout.from_assembly_graph(assembly_graph).paths``.
    """
    return LinearReadout.from_assembly_graph(assembly_graph).primary_sequence_parts()
