from typing import Dict, List, Optional, Set
import logging

from rdkit import Chem

from retromol.data.motifs import MOTIFS


class Motif:
    """Defines a chemical motif."""

    def __init__(self, name: str, smiles: str) -> None:
        """Initializes a Motif chemical motif.
        
        :param name: The name of the motif.
        :type name: str
        :param smiles: The SMILES pattern of the motif.
        :type smiles: str
        """
        self._name = name
        self._smiles = smiles
        self._mol = Chem.MolFromSmiles(smiles)

    @property
    def name(self) -> str:
        """The name of the motif."""
        return self._name
    
    @property
    def smiles(self) -> str:
        """The SMILES pattern of the motif."""
        return self._smiles
    
    @property
    def mol(self):
        """The RDKit molecule object."""
        return self._mol
    
    def is_match(self, mol: Chem.Mol) -> bool:
        """Checks if the motif matches the molecule exactly.
        
        :param mol: The molecule to check.
        :type mol: Chem.Mol
        :return: True if the motif matches the molecule, False otherwise.
        :rtype: bool
        """
        has_substruct_match = mol.HasSubstructMatch(self._mol)
        has_equal_num_atoms = mol.GetNumAtoms() == self._mol.GetNumAtoms()
        return has_substruct_match and has_equal_num_atoms
    

LOADED_MOTIFS = [Motif(motif["name"], motif["smiles"]) for motif in MOTIFS]


def get_default_motifs() -> List[Motif]:
    """Returns the default motifs for the retrosynthesis.
    
    :return: The default motifs for the retrosynthesis.
    :rtype: List[Motif]
    """
    return LOADED_MOTIFS


def match_mol_greedily(mol: Chem.Mol, logger: Optional[logging.Logger] = None) -> Optional[str]:
    """Match a molecule to a motif.
    
    :param mol: RDKit molecule
    :type mol: Chem.Mol
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: The name of the motif if a match is found, otherwise None
    :rtype: Optional[str]

    .. note:: This function uses a greedy approach to match a molecule to a motif.
    """
    motifs = get_default_motifs()
    for motif in motifs:
        if motif.is_match(mol):
            if logger: logger.debug(f"matched motif {motif.name} to molecule {Chem.MolToSmiles(mol)}")
            return motif.name
    
    # if no motif matched, return None
    if logger: logger.debug(f"no motif matched to molecule {Chem.MolToSmiles(mol)}")
    return None


def greedy_max_set_cover(
    encoding_to_mol: Dict[int, Chem.Mol], 
    nodes: List[int],
    logger: Optional[logging.Logger] = None
) -> List[int]:
    """Find biggest non-overlapping set of mol nodes in the reaction graph.

    :param encoding_to_mol: Mapping of encoding to RDKit molecule.
    :type encoding_to_mol: Dict[int, Chem.Mol]
    :param nodes: List of node encodings to consider for set cover.
    :type nodes: List[int]
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: List of selected node encodings.
    :rtype: List[int]
    """
    # create subsets of atom mappings per node.
    subsets = list()
    for node in nodes:
        mol = encoding_to_mol[node]
        tags = {atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetIsotope() != 0}
        subsets.append((node, tags))

    # sort subsets by size of atom mappings, from largest to smallest.
    sorted_subsets = sorted(subsets, key=lambda x: len(x[1]), reverse=True)
    if logger: logger.debug(f"found {len(sorted_subsets)} subsets to consider for set cover")

    # perform greedy set cover algorithm.
    selected_subsets: List[int] = []
    covered_elements: Set[int] = set()
    for node, subset in sorted_subsets:
        uncovered_elements = subset - covered_elements

        # make sure that only a subset is selected if all elements are uncovered.
        if uncovered_elements != subset:
            continue

        if uncovered_elements:
            selected_subsets.append(node)
            covered_elements.update(uncovered_elements)

    if logger: logger.debug(f"greedily selected {len(selected_subsets)} subsets for set cover")

    return selected_subsets