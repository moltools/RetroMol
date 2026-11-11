"""Data structures for representing readouts from RetroMol parsing results."""

from dataclasses import dataclass
from typing import Any, Callable, Literal, TypeVar

import networkx as nx

from retromol.chem.mol import smarts_to_mol
from retromol.model.reaction_graph import MolNode, ReactionGraph
from retromol.model.assembly_graph import AssemblyGraph
from retromol_fingerprint.fingerprint import TOKEN_LINK


ReadoutMode = Literal["leaf_identified", "first_identified"]

T = TypeVar("T")

_CARBOXYL = smarts_to_mol("[CX3](=[OX1])-[OX2H1,OX1-]")
_UNIT_ENTRY_CAP = smarts_to_mol("[#6]-[0S]-[0OH]")
_AMINO_ACID_ENTRY = smarts_to_mol("[NX3,NX4+]-[CX4]-[CX3](=[OX1])-[OX2H1,OX1-]")


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


def _backbone_entry_tags(node: MolNode) -> set[int]:
    """Incoming carbon of a capped unit or backbone N of an amino acid.

    The N must be attached to the carbon next to the fragment's carboxyl;
    unrelated amino groups (such as lysine's side-chain N) are not entries.
    Original tags let these same sites be located in the assembly's root bonds.
    """
    return {
        node.mol.GetAtomWithIdx(match[0]).GetIsotope()
        for pattern in (_UNIT_ENTRY_CAP, _AMINO_ACID_ENTRY)
        for match in node.mol.GetSubstructMatches(pattern)
        if node.mol.GetAtomWithIdx(match[0]).GetIsotope()
    }


def _remove_side_links(g: nx.Graph) -> None:
    """Omit side attachments that compete with connected backbone entries.

    A capped unit's incoming carbon and an amino acid's backbone nitrogen are
    equivalent entry sites. If an entry is connected to a preceding carboxyl
    side, an additional acyl link to a non-entry O/S/N on the same unit is a
    side attachment, whether it closes a ring or joins separate chains.
    Being a graph bridge does not make that attachment a backbone continuation.
    An extra ether link to the already connected incoming carbon likewise does
    not replace the connection from the preceding unit's carboxyl side.
    Acyl links without this competing entry evidence retain their connectivity.
    This mutates only the readout's graph copy, never the source assembly.
    """
    nodes = {key: data["molnode"] for key, data in g.nodes(data=True)
             if data["molnode"].identity.matched_rule.sequence_role != "tailoring"}
    entries = {key: _backbone_entry_tags(node) for key, node in nodes.items()}
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
            for carbon, hetero, carbon_symbol, hetero_symbol in (
                (bond.a1_tag, bond.a2_tag, bond.a1_symbol, bond.a2_symbol),
                (bond.a2_tag, bond.a1_tag, bond.a2_symbol, bond.a1_symbol),
            ):
                if bond.bond_type != "SINGLE" or carbon_symbol != "C":
                    continue
                for donor, acceptor in ((u, v), (v, u)):
                    if hetero not in g.nodes[acceptor]["tags"]:
                        continue
                    if (hetero_symbol in ("O", "S", "N")
                            and carbon in retained[donor]
                            and hetero not in entries[acceptor]
                            and connected_entries[acceptor]):
                        side_attachment = True
                    if hetero_symbol in ("O", "S") and carbon in connected_entries[donor]:
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


def _backbone_connections(g: nx.Graph) -> nx.DiGraph:
    """Directed continuations inferred solely from this assembly's fragments.

    Exit-to-entry bonds continue a chain. Bonds between two exits, two entries,
    or side sites of already connected units are crosslinks, not continuations.
    Missing port annotations retain the element/role fallback from _readout_graph;
    missing evidence must not be confused with contradictory evidence.

    Anhydride-style cuts can cap a segment in the opposite direction. A free
    retained acid anchors one end; a separate acyl starter must anchor the other.
    Reverse only the artificial-cap connections between these anchors;
    retained amide/ester directions never reverse. This resolves the direction
    of the supplied fragments without replacing, rematching or reparsing them.
    """
    entries = {key: _backbone_entry_tags(data["molnode"]) for key, data in g.nodes(data=True)}
    exits = {key: _carboxyl_exit_tags(data["molnode"]) for key, data in g.nodes(data=True)}
    retained = {key: _carboxyl_exit_tags(data["molnode"], include_added_caps=False)
                for key, data in g.nodes(data=True)}
    capped_entries = {
        key: {data["molnode"].mol.GetAtomWithIdx(match[0]).GetIsotope()
              for match in data["molnode"].mol.GetSubstructMatches(_UNIT_ENTRY_CAP)}
        for key, data in g.nodes(data=True)
    }
    connected = {key: {tag for _, _, data in g.edges(key, data=True)
                       for bond in data["bonds"] for tag in (bond.a1_tag, bond.a2_tag)} for key in g}
    terminals = {key for key in g if retained[key] and not retained[key] & connected[key]}
    sources = {key for key in g if exits[key] and not entries[key] and key not in terminals}
    strict = nx.DiGraph()
    strict.add_nodes_from(g.nodes(data=True))
    reversible = nx.DiGraph()
    reversible.add_nodes_from(g)
    links = []
    for u, v, data in g.edges(data=True):
        for bond in data["bonds"]:
            a, b = bond.a1_tag, bond.a2_tag
            if a not in g.nodes[u]["tags"]:
                a, b = b, a
            links.append((u, v, a, b))
            for donor, acceptor, out_tag, in_tag in ((u, v, a, b), (v, u, b, a)):
                if out_tag in exits[donor] and in_tag in entries[acceptor]:
                    strict.add_edge(donor, acceptor)
                    if out_tag not in retained[donor] and in_tag in capped_entries[acceptor]:
                        reversible.add_edge(donor, acceptor)

    reverse_edges = set()
    reversed_units = set()
    for terminal in terminals:
        if entries[terminal] or not exits[terminal] - retained[terminal]:
            continue
        segment = {terminal} | nx.descendants(reversible, terminal)
        # A free side-chain acid on a starter is insufficient. The opposite
        # end must actually meet a separate acyl starter at two outgoing ports.
        # Only this combination anchors both ends of a reverse-cut segment.
        if not any(
            donor in sources and acceptor in segment
            and out_tag in exits[donor] and in_tag in retained[acceptor]
            for u, v, a, b in links
            for donor, acceptor, out_tag, in_tag in ((u, v, a, b), (v, u, b, a))
        ):
            continue
        reversed_units.update(segment)
        reverse_edges.update(reversible.subgraph(segment).edges())
    directed = strict.copy()
    directed.remove_edges_from(reverse_edges)
    directed.add_edges_from((v, u) for u, v in reverse_edges)

    for u, v, a, b in links:
        if strict.has_edge(u, v) or strict.has_edge(v, u):
            continue
        a_exit, b_exit = a in exits[u], b in exits[v]
        if a_exit and b_exit:
            # Join a reverse-capped segment to its entryless acyl starter.
            # Two ordinary outgoing ports (e.g. dimerisation) stay separate.
            if u in sources and v in reversed_units:
                directed.add_edge(u, v)
            if v in sources and u in reversed_units:
                directed.add_edge(v, u)
            continue
        if a_exit or b_exit:
            donor, acceptor, tag = (u, v, b) if a_exit else (v, u, a)
            acceptor_mol = g.nodes[acceptor]["molnode"].mol
            hetero_entry = any(atom.GetIsotope() == tag and atom.GetAtomicNum() in (7, 8, 16)
                              for atom in acceptor_mol.GetAtoms())
            if (not (entries[acceptor] or exits[acceptor]) or hetero_entry) and strict.in_degree(acceptor) == 0:
                directed.add_edge(donor, acceptor)
            continue
        if a in entries[u] and not (entries[v] or exits[v]):
            directed.add_edge(v, u)
        elif b in entries[v] and not (entries[u] or exits[u]):
            directed.add_edge(u, v)
        elif not (entries[u] or exits[u] or entries[v] or exits[v]):
            # No chemical direction is available for these legacy/custom units.
            directed.add_edge(u, v)
            directed.add_edge(v, u)
    return directed


def _directed_paths(g: nx.DiGraph) -> list[list[MolNode]]:
    """Cover every unit once, following continuations rather than crosslinks.

    Start at unconnected entries; a source-free ring gets a conventional cut.
    Reachability lookahead avoids shortcuts that strand intermediate units.
    Branches become additional paths, never an undirected walk through a junction.
    Ties have a stable order; this does not assert a unique biosynthetic route.
    """
    g = g.copy()

    def node_key(key):
        node = g.nodes[key]["molnode"]
        return _node_name(node), tuple(sorted(g.nodes[key]["tags"])), node.enc

    def walk(start):
        path = [start]
        remaining = set(g) - {start}
        while remaining:
            options = [key for key in g.successors(path[-1]) if key in remaining]
            if not options:
                break
            rest = g.subgraph(remaining)
            next_node = min(options, key=lambda key: (
                -len(nx.descendants(rest, key)), rest.in_degree(key), node_key(key),
            ))
            path.append(next_node)
            remaining.remove(next_node)
        return path

    paths = []
    while g:
        starts = [key for key, degree in g.in_degree() if degree == 0]
        if not starts:
            # A ring can have an outgoing chain. Start in a source strongly
            # connected component, not at an alphabetically earlier tail unit.
            condensed = nx.condensation(g)
            starts = [key for component, degree in condensed.in_degree() if degree == 0
                      for key in condensed.nodes[component]["members"]]
        candidates = [walk(key) for key in sorted(starts, key=node_key)]
        path = min(candidates, key=lambda path: (-len(path), tuple(node_key(key) for key in path)))
        paths.append([g.nodes[key]["molnode"] for key in path])
        g.remove_nodes_from(path)
    return paths


def _assembly_paths(assembly: AssemblyGraph) -> list[list[MolNode]]:
    """One canonical partition for one supplied assembly, including saved results."""
    g = _readout_graph(assembly)
    paths = []
    for component in nx.connected_components(g):
        h = g.subgraph(component).copy()
        if any(data["molnode"].identity.matched_rule.sequence_role == "tailoring"
               for _, data in h.nodes(data=True)) or not any(
                   _backbone_entry_tags(data["molnode"]) or _carboxyl_exit_tags(data["molnode"])
                   for _, data in h.nodes(data=True)):
            # Tailoring groups and units with no chemical ports keep the existing
            # deterministic topology-only convention.
            pending = [AssemblyGraph(h)]
            while pending:
                part = pending.pop()
                path = part.longest_path()
                paths.append(path)
                used = {node.enc for node in path}
                part.g.remove_nodes_from([key for key, data in part.g.nodes(data=True)
                                          if data["molnode"].enc in used])
                pending.extend(part.connected_components())
        else:
            directed = _backbone_connections(h)
            for component in nx.weakly_connected_components(directed):
                paths.extend(_directed_paths(directed.subgraph(component)))
    return sorted(paths, key=lambda path: (-len(path), [_node_name(node) for node in path],
                                          [node.enc for node in path]))


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
    :cvar paths: Derived disjoint paths of identified MolNodes. Recomputed from
        assembly_graph on construction, including when loading old saved paths.
        Unresolved regions break paths and remain only in the source graph.
    """

    assembly_graph: AssemblyGraph
    paths: list[list[MolNode]]

    def __post_init__(self) -> None:
        """The assembly is authoritative; legacy saved paths are only a cache."""
        object.__setattr__(self, "paths", _assembly_paths(self.assembly_graph))

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
        """Read one supplied graph as disjoint N/starter-to-carboxyl paths.

        Fragment connection sites separate continuations from side attachments
        and crosslinks. Every identified unit occurs exactly once. Tailoring
        groups remain separate; gaps break paths. The graph and its matched units
        are unchanged. Cyclic chains have a conventional start, and units lacking
        chemical direction use deterministic topology-only ordering.
        """
        return cls(assembly_graph=assembly_graph, paths=[])

    def primary_sequence_parts(self) -> list[list[str]]:
        """Named backbone/tailoring paths, longest first with lexical tie-breaks.

        Each identified assembly monomer occurs once, including repeated
        identities. Unresolved and unassigned regions are omitted; assembly
        coverage measures completeness. Sublists do not encode attachment
        points, omitted crosslinks, or ring closures; retain ``assembly_graph``
        when those are needed. Untagged paths are not biological role labels.
        Backbone continuations run N/starter-to-carboxyl. Rings have a
        conventional cut; unannotated units retain deterministic ordering.
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
        return cls.from_assembly_graph(AssemblyGraph.from_dict(data["assembly_graph"]))


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
