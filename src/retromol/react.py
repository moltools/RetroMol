from typing import List, Optional, Set

from rdkit import Chem
from rdkit.Chem.rdChemReactions import ReactionFromSmarts

from retromol.data.motifs import MOTIFS
from retromol.data.rules import LINEARIZATION_RULES, SEQUENCING_RULES


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


class Reaction:
    """Defines a chemical reaction."""

    def __init__(self, name: str, reaction_smarts: str) -> None:
        """Initializes a Reaction chemical reaction.
        
        :param name: The name of the reaction.
        :type name: str
        :param smarts: The SMARTS pattern of the reaction.
        :type smarts: str
        """
        self._name = name
        self._reaction_smarts = reaction_smarts
        self._rxn = ReactionFromSmarts(reaction_smarts)

    @property
    def name(self) -> str:
        """The name of the reaction."""
        return self._name
    
    @property
    def reaction_smarts(self) -> str:
        """The SMARTS pattern of the reaction."""
        return self._reaction_smarts
    
    @property
    def rxn(self):
        """The RDKit reaction object."""
        return self._rxn
    
    def __call__(self, reactant: Chem.Mol, mask: Optional[Set[int]] = None) -> List[List[Chem.Mol]]:
        """Applies the reaction to the reactant.
        
        :param reactant: The reactant molecule.
        :type reactant: Chem.Mol
        :param mask: The atom indices allowed to be modified by the reaction.
        :type mask: Optional[Set[int]]
        :return: The products of the reaction.
        :rtype: List[List[Chem.Mol]]
        """
        # TODO
        raise NotImplementedError


def get_default_motifs() -> List[Motif]:
    """Returns the default motifs for the retrosynthesis.
    
    :return: The default motifs for the retrosynthesis.
    :rtype: List[Motif]
    """
    return [Motif(motif["name"], motif["smiles"]) for motif in MOTIFS]


def get_default_linearization_rules() -> List[Reaction]:
    """Returns the default rules for the retrosynthesis.
    
    :return: The default rules for the retrosynthesis.
    :rtype: List[Reaction]
    """
    return [Reaction(rule["name"], rule["reaction_smarts"]) for rule in LINEARIZATION_RULES]


def get_default_sequencing_rules() -> List[Reaction]:
    """Returns the default rules for the retrosynthesis.
    
    :return: The default rules for the retrosynthesis.
    :rtype: List[Reaction]
    """
    return [Reaction(rule["name"], rule["reaction_smarts"]) for rule in SEQUENCING_RULES]
