# -*- coding: utf-8 -*-

"""This module contains the API for RetroMol."""

import logging
import typing as ty

from rdkit import Chem

from retromol.retrosynthesis.chem import MolecularPattern, Molecule, ReactionRule
from retromol.retrosynthesis.graph import reaction_tree_to_digraph, reaction_tree_to_monomer_graph
from retromol.retrosynthesis.result import Result
from retromol.retrosynthesis.sequencing import parse_modular_natural_product


def parse_reaction_rules(data: ty.Dict[str, str]) -> ty.List[ReactionRule]:
    """Parse reaction rules from JSON.

    :param data: The JSON source.
    :type data: str
    :return: The reaction rules.
    :rtype: ty.List[ReactionRule]
    :raises Exception: If 'name' and/or 'pattern' is not provided for JSON
        item.
    """
    reaction_rules = []
    for item in data:
        name = item.get("name", None)
        pattern = item.get("pattern", None)

        if name is None or pattern is None:
            raise Exception("Reaction item invalid: 'name' and/or 'pattern' is not " "provided.")

        reaction_rule = ReactionRule(name, pattern)
        reaction_rules.append(reaction_rule)

    return reaction_rules


def parse_molecular_patterns(data: str) -> ty.List[MolecularPattern]:
    """Parse molecular patterns from JSON.

    :param data: The JSON source.
    :type data: str
    :return: The molecular patterns.
    :rtype: ty.List[MolecularPattern]
    :raises Exception: If 'name' and/or 'pattern' is not provided for JSON
        item.
    """
    molecular_patterns = []
    for item in data:
        name = item.get("name", None)
        pattern = item.get("pattern", None)

        if name is None or pattern is None:
            raise Exception("Monomer item invalid: 'name' and/or 'pattern' is not " "provided.")

        molecular_pattern = MolecularPattern(name, pattern)
        molecular_patterns.append(molecular_pattern)

    return molecular_patterns


def parse_mol(
    mol: Molecule,
    reactions: ty.List[ReactionRule],
    monomers: ty.List[MolecularPattern],
    neutralize: bool = True,
) -> Result:
    """Parse a molecule with RetroMol.

    :param mol: The molecule.
    :type mol: Molecule
    :param reactions: The reaction rules.
    :type reactions: ty.List[ReactionRule]
    :param monomers: The molecular patterns.
    :type monomers: ty.List[MolecularPattern]
    :param neutralize: Whether to neutralize the molecule.
    :type neutralize: bool
    :return: The result.
    :rtype: Result
    """
    logger = logging.getLogger(__name__)

    # check if all reaction names are unique
    reaction_names = [reaction.name for reaction in reactions]
    if len(reaction_names) != len(set(reaction_names)):

        # get non unique reaction names
        non_unique_reaction_names = list(set([
            reaction_name 
            for reaction_name in reaction_names 
            if reaction_names.count(reaction_name) > 1
        ]))

        msg = f"Reaction names are not unique: {non_unique_reaction_names}"
        logger.error(msg)
        raise ValueError(msg)

    name = mol.name

    if neutralize:
        mol.neutralize()

    logger.debug(f"Starting parsing for {name} ...")
    reactant, reaction_tree, reaction_mapping = mol.apply_rules(reactions)
    logger.debug(f"Applied reactions to {name}.")

    applied_reactions = list(
        set([reaction for _, reactions in reaction_tree.items() for reaction in reactions])
    )
    logger.debug(f"Applied reactions: {applied_reactions}")

    logger.debug(f"Starting monomer graph generation for {name} ...")
    reaction_tree = reaction_tree_to_digraph(reaction_tree)
    logger.debug("... Created digraph from reaction tree.")
    monomer_graph, monomer_mapping = reaction_tree_to_monomer_graph(
        mol, reaction_tree, reaction_mapping, monomers
    )
    logger.debug("... Generated monomer graph from digraph.")
    logger.debug(f"Generated monomer graph for {name}.")

    # Reformat reaction_tree and reaction_mapping into one graph.
    logger.debug(f"Reformatting reaction tree for {name} ...")
    new_reaction_tree = {}
    for parent in reaction_tree.nodes:
        children = list(reaction_tree.successors(parent))

        mol = reaction_mapping[parent]
        for atom in mol.GetAtoms():
            amn = atom.GetIsotope()
            atom.SetIsotope(0)
            atom.SetAtomMapNum(amn)
        smiles = Chem.MolToSmiles(mol)

        new_reaction_tree[parent] = {"smiles": smiles, "children": children}

    logger.debug(f"Reformatted reaction tree for {name}.")

    # Reformat monomer_graph and monomer_mapping into one graph.
    logger.debug(f"Reformatting monomer graph for {name} ...")
    new_monomer_graph = {}

    for reaction_tree_node, (monomer_graph_node, identity) in monomer_mapping.items():
        new_monomer_graph[monomer_graph_node] = {
            "reaction_tree_id": reaction_tree_node,
            "identity": identity,
            "neighbors": [],
        }

    for node in monomer_graph:
        neighbors = list(monomer_graph[node])

        if node in new_monomer_graph:
            new_monomer_graph[node]["neighbors"] = neighbors
        else:
            new_monomer_graph[node] = {
                "reaction_tree_id": None,
                "identity": None,
                "neighbors": neighbors,
            }

    logger.debug(f"Reformatted monomer graph for {name}.")

    if logger.isEnabledFor(logging.DEBUG):
        logger.debug(f" Parsing was successful for {name}.")

        for node, items in new_monomer_graph.items():
            if items["identity"] is not None:
                logger.debug(f"Found identity for {node}: {items['identity']}")

    try:
        used_reactions = [r for r in reactions if r.name in applied_reactions]
        seqs = parse_modular_natural_product(new_reaction_tree, new_monomer_graph)
    except Exception as e:
        logger.error(f" Error while parsing modular natural product: {e}")
    finally:
        if logger.isEnabledFor(logging.DEBUG):
            if not len(seqs):
                logger.debug("No modular natural product sequence found.")
            for seq in seqs:
                logger.debug(f"Modular natural product sequence: {seq}")

    # Create result object.
    result = Result(
        identifier=name,
        mol=reactant,
        success=True,
        reaction_tree=new_reaction_tree,
        applied_reactions=applied_reactions,
        monomer_graph=new_monomer_graph,
        sequences=seqs,
    )

    return result
