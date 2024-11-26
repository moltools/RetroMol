
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
    # retrieve reaction rules
    if logger: logger.debug("loading default reaction rules")
    linearization_rules = get_default_linearization_rules()
    if logger: logger.debug(f"loaded {len(linearization_rules)} linearization reaction rules")
    sequencing_rules = get_default_sequencing_rules()
    if logger: logger.debug(f"loaded {len(sequencing_rules)} sequencing reaction rules")

    # parse the molecule
    reactant = Reactant(name, smiles)

    # preprocess the molecule
    if logger: logger.debug(f"preprocessing molecule {name} with linearization rules")
    encoding_to_mol, encoding_to_rxn, reaction_graph = preprocess_mol(reactant, linearization_rules, logger)
    leaf_nodes = [parent for parent, rxns in reaction_graph.items() if not rxns]
    if logger: logger.debug(f"found {len(leaf_nodes)} leaf nodes")

    # TODO: the atom tags are also not correctly reset during/after sequencing, we need this
    # to accurately map the individual motifs to the original molecule

    # find identities for nodes
    encoding_to_identity = {}
    for encoding, mol in encoding_to_mol.items():
        encoding_to_identity[encoding] = match_mol_greedily(mol, logger)

    # try to sequence all the unidentified leaf nodes
    # these are arguably the core of the molecule consisting of a sequence of motifs
    encoding_to_motif_codes = {}
    for mol_encoding in leaf_nodes:

        # skip if node is already identifiable
        if encoding_to_identity[mol_encoding] is not None:
            continue
        
        motif_codes = sequence_mol(encoding_to_mol[mol_encoding], sequencing_rules, logger)
        if motif_codes:
            encoding_to_motif_codes[mol_encoding] = motif_codes

    # find best match for identified nodes
    matched_nodes =  [n for n, identity in encoding_to_identity.items() if identity is not None]
    identified_picked_nodes = greedy_max_set_cover(encoding_to_mol, matched_nodes)

    # gather all tags of the identified nodes
    matched_tags =  set()
    for node in identified_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        matched_tags.update(node_tags)

    # identified_picked_nodes + unidentified_picked_nodes should be the set of all nodes from the input molecule
    all_tags = set(reactant.get_tags())

    # filter out any tags from all_tags that are not in matched_tags
    unmatched_tags = all_tags - matched_tags

    # Find all nodes from the reaction tree that only contain unmatched tags.
    unmatched_nodes = []
    for node, node_mol in encoding_to_mol.items():
        tags = [atom.GetIsotope() for atom in node_mol.GetAtoms() if atom.GetIsotope() != 0]
        if all(tag in unmatched_tags for tag in tags):
            unmatched_nodes.append(node)
    
    # find best match for unidentified nodes
    unidentified_picked_nodes = greedy_max_set_cover(encoding_to_mol, unmatched_nodes)

    tags_identified_nodes = set()
    for node in identified_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_identified_nodes.update(node_tags)

    tags_unidentified_nodes = set()
    for node in unidentified_picked_nodes:
        node_tags = [atom.GetIsotope() for atom in encoding_to_mol[node].GetAtoms() if atom.GetIsotope() != 0]
        tags_unidentified_nodes.update(node_tags)

    # sanity check, all tags should be either matched or unmatched
    assert not tags_identified_nodes.intersection(tags_unidentified_nodes), f"tags {tags_identified_nodes.intersection(tags_unidentified_nodes)} are both matched and unmatched at the same time"
    assert tags_identified_nodes.union(tags_unidentified_nodes) == all_tags, f"tags {all_tags - tags_identified_nodes.union(tags_unidentified_nodes)} are missing"

    # calculate coverage score, round to 3 decimal places
    coverage_score = len(tags_identified_nodes) / len(all_tags) * 100
    coverage_score = round(coverage_score, 2)

    # TODO: generate output files
    # 1) data structure that shows the monomers that the reactant molecule is composed of
    # 2) data structure that allows reconstruction of the reactant molecule

    return coverage_score