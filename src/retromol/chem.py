from typing import List, Optional, Set
import logging

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
import numpy as np

from retromol.data.motifs import MOTIFS
from retromol.data.rules import LINEARIZATION_RULES, SEQUENCING_RULES


def mol_to_fingerprint(mol: Chem.Mol, radius: int = 2, num_bits: int = 2048):
    """Converts an RDKit molecule to a Morgan fingerprint.

    :param mol: The molecule to convert.
    :type mol: Chem.Mol
    :param radius: The radius of the fingerprint.
    :type radius: int
    :param num_bits: The number of bits in the fingerprint.
    :type num_bits: int
    :return: The Morgan fingerprint.
    """
    fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=num_bits)
    return fp_generator.GetFingerprint(mol)


def encode_mol(mol: Chem.Mol, additional_bits: Optional[List[int]] = None) -> int:
    """Encodes an RDKit molecule as an integer.

    :param mol: The molecule to encode.
    :type mol: Chem.Mol
    :param additional_bits: Additional bits to include in the encoding.
    :type additional_bits: Optional[List[int]
    :return: The encoded molecule.
    :rtype: int
    """
    fp = mol_to_fingerprint(mol)
    if additional_bits is None:
        additional_bits = []
    return hash(np.hstack([fp, additional_bits]).data.tobytes())


def neutralize_mol(mol: Chem.Mol) -> None:
    """Neutralizes the charges on an RDKit molecule.
    
    :param mol: The molecule to neutralize.
    :type mol: Chem.Mol

    .. note:: This function modifies the input molecule in place.
    """
    charge_smarts = "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"
    charge_pattern = Chem.MolFromSmarts(charge_smarts)
    at_matches = mol.GetSubstructMatches(charge_pattern)
    if len(at_matches) > 0:
        for match in at_matches:
            at_idx = match[0]  # get the atom index from the match tuple
            atom = mol.GetAtomWithIdx(at_idx)
            charge = atom.GetFormalCharge()
            h_count = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(h_count - charge)
            atom.UpdatePropertyCache()


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

        # check if reaction has 1 reactant and >=1 products
        num_reactants = self._rxn.GetNumReactantTemplates()
        num_products = self._rxn.GetNumProductTemplates()
        if num_reactants != 1 or num_products < 1:
            raise ValueError(f"reaction must have 1 reactant and >=1 products: {name}")

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
    
    def __call__(
        self, 
        reactant: Chem.Mol, 
        mask: Optional[Set[int]] = None,
        logger: Optional[logging.Logger] = None,
    ) -> List[List[Chem.Mol]]:
        """Applies the reaction to the reactant.
        
        :param reactant: The reactant molecule.
        :type reactant: Chem.Mol
        :param mask: The atom tags allowed to be modified by the reaction.
        :type mask: Optional[Set[int]]
        :param logger: The logger to use for logging.
        :type logger: Optional[logging.Logger]
        :return: The products of the reaction.
        :rtype: List[List[Chem.Mol]]
        """
        # mask is set of atom tags that are allowed to be modified by the reaction
        if mask is not None:
            # check if any tag is 0
            for atom in reactant.GetAtoms():
                if atom.GetIsotope() == 0:
                    raise ValueError(f"before using reaction {self.name}: atom has tag 0")
            tag_to_atomic_num = {atom.GetIsotope(): atom.GetAtomicNum() for atom in reactant.GetAtoms()}

            # check if any atomic num is 0:
            for _, atomic_num in tag_to_atomic_num.items():
                if atomic_num == 0:
                    raise ValueError(f"before using reaction {self.name}: atom has atomic num 0")
                
            # now we can apply the mask
            for atom in reactant.GetAtoms():
                if atom.GetIsotope() not in mask:
                    atom.SetAtomicNum(0)
        
        # apply reaction
        results = self.rxn.RunReactants([reactant])
        if logger and mask: logger.debug(f"reaction {self.name} applied with mask {mask} to reactant resulted in {len(results)} products")

        # if no results, return empty list
        if not results:
            return []

        # sanitize results
        sanitization_results = [True for _ in results]
        for result_idx, result in enumerate(results):
            for product in result:
                if mask is not None:
                    # reset atomic nums
                    for atom in product.GetAtoms():
                        if atom.GetIsotope() != 0:  # newly added atoms in reaction have isotope 0
                            original_atomic_num = tag_to_atomic_num.get(atom.GetIsotope(), None)
                            if original_atomic_num is None:
                                raise ValueError(f"no atomic num found for atom tag {atom.GetIsotope()}")
                            atom.SetAtomicNum(original_atomic_num)

                    # check if any atomic num is 0:
                    for atom in product.GetAtoms():
                        if atom.GetAtomicNum() == 0:
                            raise ValueError(f"after using reaction {self.name}: atom has atomic num 0")
    
                for atom in product.GetAtoms():
                    # check if there is aromaticity, don't touch these atoms
                    if atom.GetIsAromatic():
                        continue

                    # sanity check if atom complies with valence rules, otherwise adjust explicit Hs
                    valence_bonds = sum([bond.GetValenceContrib(atom) for bond in atom.GetBonds()])
                    valence_bonds = int(valence_bonds)
                    default_valence = Chem.GetPeriodicTable().GetDefaultValence(atom.GetAtomicNum())
                    num_hs = atom.GetNumExplicitHs()
                    if default_valence - valence_bonds < num_hs:
                        new_valence = default_valence - valence_bonds

                        # if new valence is negative, set to 0
                        if new_valence < 0:
                            msg = f"new valence is negative for atom {atom.GetSymbol()} in mol {Chem.MolToSmiles(product)} after reaction {self.name}"
                            if logger: logger.error(msg)
                            raise ValueError(msg)

                        atom.SetNumExplicitHs(new_valence)

                # sanitize product
                try:
                    Chem.SanitizeMol(product, catchErrors=True)
                except ValueError:
                    if logger: logger.warning(f"failed sanitization of product {Chem.MolToSmiles(product)} in reaction {self.name}")
                    sanitization_results[result_idx] = False
                    break

        # filter out failed sanitization results
        results = [result for result, success in zip(results, sanitization_results) if success]

        # dereplicate results, because sometimes reaction template is symmetrical and produces the same product twice
        # but in a different order; since we are not interested in the order, we can remove duplicates
        if len(results) > 1:
            encoded_results = []
            for result in results:
                encoded_result = []
                for product in result:
                    encoded_result.append(encode_mol(product))
                encoded_results.append(frozenset(encoded_result))
            indices_unique_results = np.unique(encoded_results, return_index=True)[1]
            results = [results[i] for i in indices_unique_results]

        # return deduplicated and sanitized results
        return results
    

LOADED_LINEARIZATION_RULES = [Reaction(rule["name"], rule["reaction_smarts"]) for rule in LINEARIZATION_RULES]
LOADED_SEQUENCING_RULES = [Reaction(rule["name"], rule["reaction_smarts"]) for rule in SEQUENCING_RULES]


def get_default_linearization_rules() -> List[Reaction]:
    """Returns the default rules for the retrosynthesis.
    
    :return: The default rules for the retrosynthesis.
    :rtype: List[Reaction]
    """
    return LOADED_LINEARIZATION_RULES


def get_default_sequencing_rules() -> List[Reaction]:
    """Returns the default rules for the retrosynthesis.
    
    :return: The default rules for the retrosynthesis.
    :rtype: List[Reaction]
    """
    return LOADED_SEQUENCING_RULES


def get_default_motifs() -> List[Motif]:
    """Returns the default motifs for the retrosynthesis.
    
    :return: The default motifs for the retrosynthesis.
    :rtype: List[Motif]
    """
    return LOADED_MOTIFS


def match_mol(mol: Chem.Mol, logger: Optional[logging.Logger] = None) -> Optional[str]:
    """Match a molecule to a motif.
    
    :param mol: RDKit molecule
    :type mol: Chem.Mol
    :param logger: Logger object
    :type logger: Optional[logging.Logger]
    :return: The name of the motif if a match is found, otherwise None
    :rtype: Optional[str]
    """
    motifs = get_default_motifs()
    for motif in motifs:
        if motif.is_match(mol):
            if logger: logger.debug(f"matched motif {motif.name} to molecule {Chem.MolToSmiles(mol)}")
            return motif.name
    
    # if no motif matched, return None
    if logger: logger.debug(f"no motif matched to molecule {Chem.MolToSmiles(mol)}")
    return None
