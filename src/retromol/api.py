
from typing import Optional
import logging

from retromol.chem import get_default_linearization_rules, get_default_sequencing_rules, match_mol
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

    # try to sequence all the leaf nodes
    for mol_encoding in leaf_nodes:
        motif_codes = sequence_mol(encoding_to_mol[mol_encoding], sequencing_rules, logger)

    # go over encoding_to_mol and try to identify
    for mol_encoding, mol in encoding_to_mol.items():
        match_mol(mol, logger)

    return None