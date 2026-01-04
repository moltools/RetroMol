"""Module for chemical valence rules and validation."""

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdchem import GetPeriodicTable, PeriodicTable


PERIODIC_TABLE: PeriodicTable = GetPeriodicTable()


def get_default_valence(atom_num: int) -> int:
    """
    Returns the default valence for a given atom number.

    :param atom_num: the atomic number
    :return: the default valence
    """
    return PERIODIC_TABLE.GetDefaultValence(atom_num)


def correct_hydrogens(mol: Mol) -> None:
    """
    Correct explicit hydrogens on atoms based on valence rules.
    
    :param mol: the RDKit molecule
    .. note:: this function modifies the molecule in place
    """
    for atom in mol.GetAtoms():
        # Skip aromatic and charged atoms
        if atom.GetIsAromatic() or atom.GetFormalCharge() != 0:
            continue

        # Skip phosphorus and sulfur (can have expanded valence)
        if atom.GetAtomicNum() in {15, 16}:
            continue

        # Check if atom complies with valence rules, otherwise adjust explicit hydrogens
        valence_bonds = int(sum([bond.GetValenceContrib(atom) for bond in atom.GetBonds()]))
        default_valence = get_default_valence(atom.GetAtomicNum())
        num_hs = atom.GetNumExplicitHs()

        if default_valence - valence_bonds < num_hs:
            new_valence = default_valence - valence_bonds

            if new_valence < 0:
                raise ValueError("new valence for atom is negative")
            
            atom.SetNumExplicitHs(new_valence)
