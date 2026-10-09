"""Data structures for representing readouts from RetroMol parsing results."""

from dataclasses import dataclass
from typing import Callable, Literal, TypeVar

from retromol.model.reaction_graph import MolNode, ReactionGraph
from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.rules import MatchingRule
from retromol.chem.mol import encode_mol
from retromol.chem.tagging import get_tags_mol
from retromol_fingerprint.fingerprint import TOKEN_LINK


ReadoutMode = Literal["leaf_identified", "first_identified"]

T = TypeVar("T")


def merge_named_paths(paths: list[list[T]], key: Callable[[T], str], link_item: T) -> list[T]:
    """
    Merge candidate primary sequence paths into a single sequence, joined by
    `link_item`, so a compound/BGC with more than one disconnected path (e.g. a
    branched/cyclic assembly, or a sugar attached via a glycosidic bond that
    AssemblyGraph never treats as a "connection") still gets exactly one primary
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
    :cvar paths: list of paths (each path is a list of MolNodes) extracted from the assembly graph, representing linear sequences of monomers.
    """

    assembly_graph: AssemblyGraph
    paths: list[list[MolNode]]

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
        Create a LinearReadout from a Result object.

        :param root_enc: Encoding of the root molecule.
        :param reaction_graph: ReactionGraph object.
        :param identified_only: Whether to only include identified monomers.
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

        # Break bonds between monomers that are not backbone-related bonds (i.e., keep C-C and C-N bonds only)
        f = a.filtered_by_root_bond_elements(allow_pairs={frozenset(("C", "C")), frozenset(("C", "N"))}, drop_isolated=False)

        # Get individual connected components from the assembly graph and extract longest paths (allow to visit each edge only once)
        hs = f.connected_components()

        paths: list[list[MolNode]] = []
        for h in hs:
            path = h.longest_path()
            paths.append(path)

        return cls(assembly_graph=a, paths=paths)

    def primary_sequence(self) -> list[str]:
        """
        The single primary sequence for this readout: every path in `self.paths`
        -- including single-node ones, e.g. a lone tailoring event like
        glycosylation or methylation that AssemblyGraph never connects to the
        main chain (it only keeps C-C/C-N bonds as "connections", so e.g. a sugar
        attached via a glycosidic C-O-C linkage always ends up disconnected) --
        merged into one sequence via `merge_named_paths` (longest first, ties
        broken lexicographically, joined by TOKEN_LINK). Nothing found in the
        assembly graph is dropped; a caller that wants just one biosynthetic
        chain, without tailoring events mixed in, can split the result back
        apart on TOKEN_LINK.

        An unidentified node is named "X", the convention used everywhere else in
        RetroMol.

        :return: the single merged primary sequence
        """
        named_paths = [
            [node.identity.matched_rule.name if node.is_identified else "X" for node in path]
            for path in self.paths
        ]
        return merge_named_paths(named_paths, key=lambda n: n, link_item=TOKEN_LINK)

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
