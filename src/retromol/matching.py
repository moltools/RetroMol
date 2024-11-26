from typing import List, Optional, Set
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
