
from typing import Optional
import logging

from retromol.chem import get_default_linearization_rules, get_default_sequencing_rules
from retromol.matching import match_mol_greedily
from retromol.react import preprocess_mol, sequence_mol


def run_retromol(name: str, smiles: str, logger: Optional[logging.Logger] = None) -> None:
    """Parse a SMILES string.
    
    :param name: Name of the molecule
    :type name: str
    :param smiles: SMILES string
    :type smiles: str
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: None
    """
    # retrieve reaction rules
    if logger: logger.debug("loading default reaction rules")
    linearization_rules = get_default_linearization_rules()
    if logger: logger.debug(f"loaded {len(linearization_rules)} linearization reaction rules")
    sequencing_rules = get_default_sequencing_rules()
    if logger: logger.debug(f"loaded {len(sequencing_rules)} sequencing reaction rules")

    # preprocess the molecule
    if logger: logger.debug(f"preprocessing molecule {name} with linearization rules")
    encoding_to_mol, encoding_to_rxn, reaction_graph = preprocess_mol(smiles, linearization_rules, logger)
    leaf_nodes = [parent for parent, rxns in reaction_graph.items() if not rxns]
    if logger: logger.debug(f"found {len(leaf_nodes)} leaf nodes")

    # find identities for nodes
    encoding_to_identity = {}
    for encoding, mol in encoding_to_mol.items():
        encoding_to_identity[encoding] = [match_mol_greedily(mol, logger)]

    # try to sequence all the unidentified leaf nodes
    # these are arguably the core of the molecule consisting of a sequence of motifs
    for mol_encoding in leaf_nodes:

        # skip if node is already identified
        if encoding_to_identity[mol_encoding] is not None:
            continue
        
        motif_codes = sequence_mol(encoding_to_mol[mol_encoding], sequencing_rules, logger)
        # TODO: find best set cover

    return None