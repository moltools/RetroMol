"""This module contains functions for graph operations with networkx."""

from collections.abc import Generator
from copy import deepcopy
from typing import Any

import networkx as nx
from networkx import Graph

import retromol.chem as chem


def mol_to_graph(mol: chem.Mol, use_tags: bool = False) -> "Graph[int | str]":
    """
    Convert RDKit molecule to networkx graph based on atom indices.

    :param mol: RDKit molecule
    :return: networkx graph
    .. note:: Nodes are atom indices (or atom tags if `use_tags` is True)
    """
    smiles = chem.mol_to_smiles(mol)
    graph: Graph[int | str] = Graph(
        smiles=smiles,
        smiles_no_tags=chem.mol_to_smiles(deepcopy(mol), remove_tags=True),
    )

    # If use_tags is True, we will use atom isotopes as tags for nodes
    for atom in mol.GetAtoms():
        if use_tags:
            atom_tag: int = atom.GetIsotope()
            if atom_tag == 0:
                atom_idx: int = atom.GetIdx()
                atom_tag = -1 * atom_idx
            graph.add_node(atom_tag)
        else:
            graph.add_node(atom.GetIdx())

    # Add edges between atoms based on bonds
    # If use_tags is True, we will use atom isotopes as tags for edges
    for bond in mol.GetBonds():
        if use_tags:
            begin_atom_idx = bond.GetBeginAtomIdx()
            begin_atom_tag = bond.GetBeginAtom().GetIsotope()
            if begin_atom_tag == 0:
                begin_atom_tag = -1 * begin_atom_idx
            end_atom_idx = bond.GetEndAtomIdx()
            end_atom_tag = bond.GetEndAtom().GetIsotope()
            if end_atom_tag == 0:
                end_atom_tag = -1 * end_atom_idx
            graph.add_edge(begin_atom_tag, end_atom_tag)
        else:
            graph.add_edge(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())

    return graph


def merge_nodes(
    graph: "Graph[int | str]",
    merged_node_id: int | str,
    nodes: list[int],
    props: dict[str, Any] | None = None,
) -> None:
    """
    Merge `nodes` into a single node `merged_node_id`.

    - Internal edges among `nodes` are removed.
    - Edges from `nodes` to outside nodes are rewired to `merged_node_id`.
    - If `merged_node_id` already exists and is in `nodes`, it is kept;
      otherwise a new node is created.
    - If `props` is provided, these attributes are set on the merged node
      (overwriting any existing attributes with the same keys).

    :param graph: the graph to modify
    :param merged_node_id: ID to keep/create for the merged node
    :param nodes: node IDs to merge (must exist in `graph`)
    :param props: attributes to assign to the merged node after merge
    .. note:: this function modifies `graph` in place and returns None
    """
    if not nodes:
        return

    # Validate nodes exist
    missing = [n for n in nodes if n not in graph]
    if missing:
        raise ValueError(f"Cannot merge {missing} - they are not in the graph.")

    nodes_set: set[int | str] = set(nodes)

    # Collect outside neighbors
    outside_neighbors: set[int | str] = set()
    for n in nodes:
        for nbr in graph.neighbors(n):
            if nbr not in nodes_set:
                outside_neighbors.add(nbr)

    # Merge paths
    if merged_node_id in nodes_set:
        # Case A: keep existing merged_node_id; remove others
        to_remove = nodes_set - {merged_node_id}
        graph.remove_nodes_from(to_remove)

        # Reattach edges to outside neighbors
        for nbr in outside_neighbors:
            graph.add_edge(merged_node_id, nbr)

    else:
        # Case B: create a brand-new merged_node_id
        if merged_node_id in graph:
            raise ValueError(
                f"Cannot create merged node {merged_node_id} because it already exists "
                "in the graph (and was not part of the nodes list)."
            )

        # Remove all original nodes, then add merged node
        graph.remove_nodes_from(nodes)
        graph.add_node(merged_node_id)

        # Reattach edges to outside neighbors
        for nbr in outside_neighbors:
            graph.add_edge(merged_node_id, nbr)

    # Apply/override attributes on the merged node
    if props:
        graph.nodes[merged_node_id].update(props)


def is_linear_graph(g: "Graph[int | str]") -> bool:
    """
    Return True if `g` is a single path (a “string” of nodes):
      - connected,
      - acyclic,
      - exactly two endpoints of degree 1 (or a single node),
      - all other nodes of degree 2.

    :param g: the graph to check
    :return: True if `g` is a linear graph, False otherwise
    """
    n = g.number_of_nodes()
    # Empty graph -> not a chain
    if n == 0:
        return False

    # Single node -> trivially linear
    if n == 1:
        return True

    # Must be connected
    if not nx.is_connected(g):
        return False

    # Tree check: exactly n-1 edges for acyclic connected graph
    if g.number_of_edges() != n - 1:
        return False

    # Count degrees
    degrees: dict[int | str, int] = dict(g.degree())
    degs: list[int] = list(degrees.values())
    num_deg1 = sum(1 for d in degs if d == 1)
    num_deg2 = sum(1 for d in degs if d == 2)

    # Exactly two endpoints (degree 1) and the rest degree 2
    return (num_deg1 == 2) and (num_deg2 == n - 2)


def get_linear_path(g: "Graph[int | str]") -> list[int | str] | None:
    """
    If `g` is a linear graph (as per is_linear_graph), return the list
    of nodes in path order; otherwise return None.

    :param g: the graph to extract the path from
    :return: list of nodes in path order, or None if not linear
    """
    # Single node case
    if g.number_of_nodes() == 1:
        return list(g.nodes())

    if is_linear_graph(g):
        # Find the two endpoints (degree == 1)
        endpoints = [node for node, deg in g.degree() if deg == 1]
        start, end = endpoints

        # The graph is a tree, so shortest_path is the unique path
        return nx.shortest_path(g, source=start, target=end)


def iter_kmers(G: "Graph[int | str]", k: int) -> Generator[tuple[int | str, ...], None, None]:
    """
    Generate all length-k node walks (k-mers) from graph G, returning node identifiers.

    Each yielded k-mer is a tuple of node identifiers (length == k).

    :param G: the graph to traverse
    :param k: the length of the k-mers to generate
    :return: a generator yielding k-mers as tuples of node identifiers
    .. note:: enumerates walks (nodes may repeat), not simple paths
    .. note:: for k == 1, yields one k-mer per node
    """
    if k < 1:
        raise ValueError("k must be >= 1")

    # k == 1: one k-mer per node
    if k == 1:
        for n in G.nodes(data=False):
            yield (n,)
        return

    # Build paths as lists of NodeT; convert to tuple only when yielding.
    stack: list[tuple[int | str, list[int | str]]] = [(start, [start]) for start in G.nodes(data=False)]

    while stack:
        node, path = stack.pop()
        if len(path) == k:
            yield tuple(path)
            continue
        for nbr in G.neighbors(node):
            # Append neighbor to the current path (nodes may repeat)
            stack.append((nbr, path + [nbr]))
