"""Module for RDKit molecule utilities."""

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.MolStandardize.rdMolStandardize import TautomerEnumerator, Uncharger
from rdkit.Chem.rdmolfiles import MolFromSmiles, MolFromSmarts, MolToSmiles
from rdkit.Chem.rdmolops import GetMolFrags, RemoveStereochemistry, SanitizeMol
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.rdmolops import (
    AssignAtomChiralTagsFromStructure,
    AssignStereochemistry,
    SetBondStereoFromDirections,
)

from retromol.chem.tagging import remove_tags
from retromol.chem.valence import correct_hydrogens
from retromol.utils.hashing import sha256_hex


def sanitize_mol(mol: Mol, fix_hydrogens: bool = False) -> bool:
    """
    Sanitizes an RDKit molecule in place, returning success status.

    :param mol: the molecule to sanitize
    :param fix_hydrogens: whether to correct hydrogen counts before sanitization
    :return: True if sanitization was successful, False otherwise
    .. note:: this function mutates the input molecule in place
    """
    try:
        if fix_hydrogens:
            correct_hydrogens(mol)
        SanitizeMol(mol)
        return True
    except Exception:
        return False
    

def reassign_stereochemistry(mol: Mol) -> Mol:
    """
    Reassign stereochemistry of a molecule without changing its identity.

    :param mol: input molecule
    :return: molecule with reassigned stereochemistry
    """
    mm = Mol(mol)
    
    SetBondStereoFromDirections(mm)

    if mm.GetNumConformers() > 0:
        # When conformers are present, reassign chiral tags from structure
        AssignAtomChiralTagsFromStructure(mm, replaceExistingTags=True)

    AssignStereochemistry(mm, cleanIt=True, force=True, flagPossibleStereoCenters=True)

    return mm


def smiles_to_mol(smiles: str) -> Mol:
    """
    Converts a SMILES string to an RDKit molecule.

    :param smiles: the SMILES string to convert
    :return: the RDKit molecule
    :raises ValueError: if the SMILES is invalid
    """
    mol = MolFromSmiles(smiles)

    if mol is None:
        raise ValueError(f"invalid SMILES: {smiles}")
    
    return mol


def smarts_to_mol(smarts: str) -> Mol:
    """
    Converts a SMARTS string to an RDKit molecule.

    :param smarts: the SMARTS string to convert
    :return: the RDKit molecule
    :raises ValueError: if the SMARTS pattern is invalid
    """
    mol: Mol | None = MolFromSmarts(smarts)

    if mol is None:
        raise ValueError(f"invalid SMARTS: {smarts}")

    return mol


def get_fragments(mol: Mol) -> tuple[Mol, ...]:
    """
    Returns the fragments of a molecule.

    :param mol: the molecule to analyze
    :return: a tuple of fragment molecules
    """
    frags: tuple[Mol, ...] = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    return frags


def count_fragments(mol: Mol) -> int:
    """
    Counts the number of fragments in a molecule.

    :param mol: the molecule to analyze
    :return: the number of fragments
    """
    return len(get_fragments(mol))


def largest_fragment(mol: Mol) -> Mol:
    """
    Return the largest fragment of a molecule (by atom count).

    :param mol: input molecule
    :return: largest fragment molecule
    """
    frags = get_fragments(mol)
    return max(frags, key=lambda m: m.GetNumAtoms()) if frags else mol


def standardize_from_smiles(
    smi: str,
    keep_stereo: bool = False,
    neutralize: bool = True,
    tautomer_canon: bool = True,
) -> Mol | None:
    """
    Standardize a molecule from its SMILES representation.

    :param smi: input SMILES string
    :param keep_stereo: whether to retain stereochemistry
    :param neutralize: whether to neutralize charges
    :param tautomer_canon: whether to canonicalize tautomers
    :return: standardized molecule or None if input SMILES is invalid
    """
    mol = smiles_to_mol(smi)
    mol = largest_fragment(mol)

    if neutralize:
        mol = Uncharger().uncharge(mol)

    if tautomer_canon:
        mol = TautomerEnumerator().Canonicalize(mol)

    sanitize_mol(mol, fix_hydrogens=False)

    if not keep_stereo:
        RemoveStereochemistry(mol)

    return mol


def mol_to_smiles(
    mol: Mol,
    include_tags: bool = False,
    isomeric: bool = True,
    canonical: bool = True
) -> str:
    """
    Converts an RDKit molecule to a SMILES string.

    :param mol: the molecule to convert
    :param include_tags: whether to include tagging information in the SMILES
    :param isomeric: whether to include isomeric information in the SMILES
    :param canonical: whether to generate a canonical SMILES
    :return: the SMILES string
    """
    if not include_tags:
        mol = remove_tags(mol, in_place=False)

    return MolToSmiles(mol, isomericSmiles=isomeric, canonical=canonical)


def mol_to_inchikey(mol: Mol) -> str:
    """
    Converts an RDKit molecule to an InChIKey.

    :param mol: the molecule to convert
    :return: the InChIKey
    """
    return MolToInchiKey(mol)


def encode_mol(mol: Mol) -> str:
    """
    Encodes an RDKit molecule as a canonical isomeric SMILES string.

    :param mol: the molecule to encode
    :return: the encoded molecule
    """
    smiles = mol_to_smiles(mol, include_tags=True, isomeric=True, canonical=True)

    return sha256_hex(smiles)
