"""Module for RDKit atom tagging utilities."""

from rdkit.Chem.rdchem import Mol


def tag_mol(mol: Mol) -> None:
    """
    Tags the atoms in an RDKit molecule with unique isotope-based tags.

    :param mol: the molecule to tag
    .. note:: this function modifies the input molecule in place
    """
    for i, atom in enumerate(mol.GetAtoms(), start=1):
        atom.SetIsotope(i)


def remove_tags(mol: Mol, in_place: bool = False) -> Mol:
    """
    Removes atom tags from an RDKit molecule.

    :param mol: the molecule to process
    :param in_place: whether to modify the input molecule in place
    :return: the molecule without atom tags
    """
    if not in_place:
        mol = Mol(mol)

    for atom in mol.GetAtoms():
        atom.SetIsotope(0)

    return mol


def get_tags_mol(mol: Mol) -> set[int]:
    """
    Get the atom tags from a molecule.

    :param mol: the molecule
    :return: unordered set of atom tags
    """
    tags: set[int] = set()
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:
            tags.add(atom.GetIsotope())

    return tags


def get_tags_mols(mols: list[Mol]) -> set[int]:
    """
    Get the atom tags from a list of molecules.

    :param mols: the list of molecules
    :return: unordered set of atom tags
    """
    tags: set[int] = set()
    for mol in mols:
        tags.update(get_tags_mol(mol))

    return tags


def all_atoms_have_unique_tags(mol: Mol) -> bool:
    """
    Check if all atoms in a molecule have unique tags.

    :param mol: the molecule
    :return: True if all atoms have unique tags, False otherwise
    """
    num_atoms = mol.GetNumAtoms()
    curr_tags = get_tags_mol(mol)

    return len(curr_tags) == num_atoms
