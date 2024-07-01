# -*- coding: utf-8 -*-

"""This module contains functions for converting reaction trees to graphs."""

import logging
import typing as ty

import networkx as nx

from retromol.retrosynthesis.chem import (
    MolecularPattern,
    Molecule,
    MonomerGraphMapping,
    ReactionTreeMapping,
    Tree,
    greedy_max_set_cover,
    identify_mol,
)


def reaction_tree_to_digraph(tree: Tree) -> nx.DiGraph:
    """Convert reaction tree to directed graph.

    :param tree: Reaction tree.
    :type tree: Tree
    :return: Directed graph.
    :rtype: nx.DiGraph
    """
    digraph = nx.DiGraph()

    for parent, reactions in tree.items():
        digraph.add_node(parent)

        for reaction, children in reactions.items():
            for child in children:

                for product in child:
                    digraph.add_node(product)
                    digraph.add_edge(parent, product, reaction=reaction)

    return digraph


def reaction_tree_to_monomer_graph(
    mol: Molecule,
    tree: nx.DiGraph,
    mapping: ReactionTreeMapping,
    monomers: ty.List[MolecularPattern],
) -> ty.Tuple[nx.Graph, MonomerGraphMapping]:
    """Convert reaction tree to monomer graph.

    :param mol: Molecule.
    :type mol: Molecule
    :param tree: Reaction tree.
    :type tree: nx.DiGraph
    :param mapping: Reaction tree mapping.
    :type mapping: ReactionTreeMapping
    :param monomers: Monomers.
    :type monomers: ty.List[MolecularPattern]
    :return: Monomer graph and monomer graph mapping.
    :rtype: ty.Tuple[nx.Graph, MonomerGraphMapping]
    """
    logger = logging.getLogger(__name__)

    logger.debug("Identifying monomers in reaction tree ...")
    identified = []
    for node in tree.nodes:
        if identity := identify_mol(mapping[node], monomers):
            identified.append((node, identity))
    logger.debug(f"Identified {len(identified)} monomers in reaction tree.")

    logger.debug("Applying greedy maximum set cover ...")
    subgraphs = greedy_max_set_cover(mol.compiled, identified, mapping)
    logger.debug(f"Applied greedy maximum set cover and found {len(subgraphs)} subgraphs.")

    # All atoms in the substrate have a valid atom map number as isotope.
    logger.debug("Creating monomer graph ...")
    monomer_graph = nx.Graph()
    for atom in mol.compiled.GetAtoms():
        monomer_graph.add_node(atom.GetIsotope())
    for bond in mol.compiled.GetBonds():
        monomer_graph.add_edge(bond.GetBeginAtom().GetIsotope(), bond.GetEndAtom().GetIsotope())
    logger.debug("... converted compiled input to networkx graph.")

    monomer_graph_mapping = dict()
    for i, subgraph in enumerate(subgraphs):
        logger.debug(f"Creating monomer graph node for subgraph {i + 1}/{len(subgraphs)} ...")

        amns = [
            atom.GetIsotope() for atom in mapping[subgraph[0]].GetAtoms() if atom.GetIsotope() > 0
        ]
        logger.debug(f"Identified AMNs in subgraph: {amns} ...")

        # Merge subgraph nodes in monomer graph.
        monomer_graph.add_node(amns[0])
        for amn in amns[1:]:
            monomer_graph = nx.contracted_nodes(monomer_graph, amns[0], amn, self_loops=False)
        logger.debug("... merged nodes with identifier AMNs in monomer graph.")

        monomer_graph_mapping[subgraph[0]] = (amns[0], subgraph[1])
    logger.debug("... created monomer graph.")
    logger.debug(f"Created monomer graph with {monomer_graph.number_of_nodes()} nodes.")

    return monomer_graph, monomer_graph_mapping
