
from typing import Optional
import logging

from rdkit import Chem

from retromol.chem import Reactant, get_default_linearization_rules, get_default_sequencing_rules
from retromol.matching import match_mol_greedily, greedy_max_set_cover
from retromol.react import preprocess_mol, sequence_mol


def run_retromol(name: str, smiles: str, logger: Optional[logging.Logger] = None) -> float:
    """Parse a SMILES string.
    
    :param name: Name of the molecule
    :type name: str
    :param smiles: SMILES string
    :type smiles: str
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: coverage score
    :rtype: float
    """
    # parse the molecule
    reactant = Reactant(name, smiles)
    if logger: logger.debug(f"parsing molecule {reactant.name} with SMILES {Chem.MolToSmiles(reactant.mol)}")

    # retrieve reaction rules
    if logger: logger.debug("loading default reaction rules")
    linearization_rules = get_default_linearization_rules()
    if logger: logger.debug(f"loaded {len(linearization_rules)} linearization reaction rules")
    sequencing_rules = get_default_sequencing_rules()
    if logger: logger.debug(f"loaded {len(sequencing_rules)} sequencing reaction rules")

    # preprocess the molecule
    if logger: logger.debug(f"preprocessing molecule {name} with linearization rules")
    encoding_to_mol, encoding_to_rxn, reaction_graph = preprocess_mol(reactant, linearization_rules, logger)
    leaf_nodes = [parent for parent, rxns in reaction_graph.items() if not rxns]
    if logger: logger.debug(f"found {len(leaf_nodes)} leaf nodes")

    # find identities for nodes
    encoding_to_identity = {}
    for encoding, mol in encoding_to_mol.items():
        mol_identity = match_mol_greedily(mol, logger)
        if mol_identity is not None:
            encoding_to_identity[encoding] = mol_identity
            continue

    # find best match for identified nodes
    matched_nodes =  [n for n, identity in encoding_to_identity.items() if identity is not None]
    identified_picked_nodes = greedy_max_set_cover(encoding_to_mol, matched_nodes)
    if logger: logger.debug(f"identified picked nodes: {identified_picked_nodes}")

    # gather all tags of the identified nodes
    matched_tags =  set()
    for node in identified_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        if logger: logger.debug(f"{node}, {Chem.MolToSmiles(encoding_to_mol[node])}, {node_tags}")
        matched_tags.update(node_tags)
        continue

    # identified_picked_nodes + unidentified_picked_nodes should be the set of all nodes from the input molecule
    all_tags = set(reactant.get_tags())

    # filter out any tags from all_tags that are not in matched_tags
    unmatched_tags = all_tags - matched_tags

    # find all nodes from the reaction tree that only contain unmatched tags, and could be sequenced
    unmatched_sequenceable_nodes = []
    encoding_to_motif_codes = {}
    for node, node_mol in encoding_to_mol.items():

        # only consider leaf nodes, as these should be the sequenceable nodes
        if node not in leaf_nodes:
            continue

        motif_codes = sequence_mol(node_mol, sequencing_rules, logger)
        encoding_to_motif_codes[node] = motif_codes
        if len(motif_codes) != 0:
            tags = [atom.GetIsotope() for atom in node_mol.GetAtoms() if atom.GetIsotope() != 0]
            if all(tag in unmatched_tags for tag in tags):
                unmatched_sequenceable_nodes.append(node)
    
    # find best match for unidentified but sequenceable nodes
    unidentified_sequenceable_picked_nodes = greedy_max_set_cover(encoding_to_mol, unmatched_sequenceable_nodes)

    # gather all tags from unidentified but sequenceable nodes
    tags_unidentified_nodes = set()
    for node in unidentified_sequenceable_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_unidentified_nodes.update(node_tags)
    unmatched_tags = unmatched_tags - tags_unidentified_nodes

    # find all nodes from the reaction that only contain new unmatched tags, and cannot be sequenced
    unmatched_unsequenceable_nodes = []
    for node, node_mol in encoding_to_mol.items():
        tags = [atom.GetIsotope() for atom in node_mol.GetAtoms() if atom.GetIsotope() != 0]
        if all(tag in unmatched_tags for tag in tags):
            if node not in unmatched_sequenceable_nodes:
                unmatched_unsequenceable_nodes.append(node)

    # find best match for unidentified and unsequenceable nodes
    if unmatched_unsequenceable_nodes:
        unidentified_unsequenceable_picked_nodes = greedy_max_set_cover(encoding_to_mol, unmatched_unsequenceable_nodes)
    else:
        unidentified_unsequenceable_picked_nodes = []

    # collate all tags in identified picked nodes and separately in unidentified picked nodes
    tags_identified_nodes = set()
    for node in identified_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_identified_nodes.update(node_tags)
    if logger:
        logger.debug(f"picked identified nodes: {identified_picked_nodes}")
        logger.debug(f"tags identified nodes: {tags_identified_nodes}")
        for node in identified_picked_nodes:
            logger.debug(f"{node} - {Chem.MolToSmiles(encoding_to_mol[node])} - {encoding_to_identity[node]}")
    
    tags_unidentified_sequenceable_nodes = set()
    for node in unidentified_sequenceable_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_unidentified_sequenceable_nodes.update(node_tags)
    if logger:
        logger.debug(f"picked unidentified sequenceable nodes: {unidentified_sequenceable_picked_nodes}")
        logger.debug(f"tags unidentified sequenceable nodes: {tags_unidentified_sequenceable_nodes}")
        for node in unidentified_sequenceable_picked_nodes:
            logger.debug(f"{node} - {Chem.MolToSmiles(encoding_to_mol[node])}")

    tags_unidentified_unsequenceable_nodes = set()
    for node in unidentified_unsequenceable_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_unidentified_unsequenceable_nodes.update(node_tags)
    if logger:
        logger.debug(f"picked unidentified unsequenceable nodes: {unidentified_unsequenceable_picked_nodes}")
        logger.debug(f"tags unidentified unsequenceable nodes: {tags_unidentified_unsequenceable_nodes}")
        for node in unidentified_unsequenceable_picked_nodes:
            logger.debug(f"{node} - {Chem.MolToSmiles(encoding_to_mol[node])}")

    # sanity checks based on tags, none of the tags sets can overlap
    assert len(tags_identified_nodes & tags_unidentified_sequenceable_nodes) == 0, "identified and unidentified sequenceable tags overlap"
    assert len(tags_identified_nodes & tags_unidentified_unsequenceable_nodes) == 0, "identified and unidentified unsequenceable tags overlap"

    # sanity checks based on tags, all sets combined should have same tags as all_tags
    assert tags_identified_nodes | tags_unidentified_sequenceable_nodes | tags_unidentified_unsequenceable_nodes == all_tags, "picked nodes have different tags than all tags"

    # find best motif coverage for the unidentified sequenceable nodes
    tags_identified_motifs = set()
    for node in unidentified_sequenceable_picked_nodes:
        motif_codes = encoding_to_motif_codes[node]
        tags_identified_motifs_motif_code = set()
        for motif_code in motif_codes:
            for motif in motif_code:
                motif_identity = match_mol_greedily(motif, logger)
                if motif_identity is not None:
                    motif_tags = [atom.GetIsotope() for atom in motif.GetAtoms() if atom.GetIsotope() != 0]
                    tags_identified_motifs_motif_code.update(motif_tags)
        if len(tags_identified_motifs_motif_code) > len(tags_identified_motifs):
            tags_identified_motifs = tags_identified_motifs_motif_code

    # calculate coverage score, round to 3 decimal places
    coverage_score = (len(tags_identified_nodes) + len(tags_identified_motifs)) / len(all_tags) * 100
    coverage_score = round(coverage_score, 2)

    # TODO: generate output files
    # 1) data structure that shows the monomers that the reactant molecule is composed of
    # 2) data structure that allows reconstruction of the reactant molecule

    return coverage_score