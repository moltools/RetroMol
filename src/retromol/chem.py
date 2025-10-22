"""This module contains functions for chemical operations using RDKit."""

import numpy as np
from numpy.typing import NDArray
from rdkit import Chem, RDLogger
from rdkit.Chem.inchi import MolToInchiKey
from rdkit.Chem.MolStandardize.rdMolStandardize import LargestFragmentChooser, TautomerEnumerator, Uncharger
from rdkit.Chem.rdchem import Atom, BondStereo, BondType, GetPeriodicTable, Mol, PeriodicTable
from rdkit.Chem.rdChemReactions import ChemicalReaction, ReactionFromSmarts
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.Chem.rdMolDescriptors import GetMorganFingerprintAsBitVect
from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromSmiles, MolToSmarts, MolToSmiles
from rdkit.Chem.rdmolops import (
    AssignStereochemistry,
    GetMolFrags,
    RemoveStereochemistry,
    SanitizeMol,
    SetBondStereoFromDirections,
)
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, TanimotoSimilarity

RDLogger.DisableLog("rdApp.*")

# Type aliases for imported functions/classes
__all__ = ["Atom", "Mol", "ChemicalReaction", "ExplicitBitVect"]


def get_periodic_table() -> PeriodicTable:
    """
    Returns the RDKit periodic table.

    :return: the RDKit periodic table
    """
    return GetPeriodicTable()


def sanitize_mol(mol: Mol) -> None:
    """
    Sanitizes an RDKit molecule.

    :param mol: the molecule to sanitize
    .. note:: this function modifies the input molecule in place
    """
    SanitizeMol(mol)


def smiles_to_mol(smiles: str, retain_largest_fragment: bool = False) -> Mol:
    """
    Converts a SMILES string to an RDKit molecule.

    :param smiles: the SMILES string to convert
    :param retain_largest_fragment: if True, retains only the largest fragment of the molecule
    :return: the RDKit molecule
    :raises ValueError: if the SMILES is invalid
    """
    mol: Mol | None = MolFromSmiles(smiles)

    match mol:
        case None:
            raise ValueError(f"invalid SMILES: {smiles}")
        case m:
            if retain_largest_fragment:
                chooser: LargestFragmentChooser = LargestFragmentChooser()
                m: Mol = chooser.choose(m)
                SanitizeMol(m)
            return m


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
    sanitize_mol(mol)
    if not keep_stereo:
        RemoveStereochemistry(mol)
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


def smarts_to_reaction(smarts: str, use_smiles: bool = False) -> ChemicalReaction:
    """
    Converts a SMARTS string to an RDKit reaction.

    :param smarts: the SMARTS string to convert
    :param use_smiles: whether to interpret the SMARTS as SMILES
    :return: the RDKit reaction
    :raises ValueError: if the SMARTS pattern is invalid
    """
    rxn: ChemicalReaction | None = ReactionFromSmarts(smarts, useSmiles=use_smiles)

    if rxn is None:
        raise ValueError(f"invalid reaction SMARTS: {smarts}")

    return rxn


def mol_to_smiles(mol: Mol, remove_tags: bool = False, isomeric: bool = True, canonical: bool = True) -> str:
    """
    Converts an RDKit molecule to a SMILES string.

    :param mol: the molecule to convert
    :param remove_tags: whether to remove atom tags (isotopes) from the SMILES
    :param isomeric: whether to include isomeric information in the SMILES
    :param canonical: whether to generate a canonical SMILES
    :return: the SMILES string
    """
    if remove_tags:
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)

    return MolToSmiles(mol, isomericSmiles=isomeric, canonical=canonical)


def mol_to_smarts(mol: Mol) -> str:
    """
    Converts an RDKit molecule to a SMARTS string.

    :param mol: the molecule to convert
    :return: the SMARTS string
    """
    return MolToSmarts(mol)


def mol_to_fpr(mol: Mol, rad: int = 2, nbs: int = 2048) -> NDArray[np.int8]:
    """
    Converts an RDKit molecule to a Morgan fingerprint.

    :param mol: the molecule to convert
    :param rad: the radius of the Morgan fingerprint
    :param nbs: the number of bits in the fingerprint
    :return: the Morgan fingerprint
    """
    gen = GetMorganGenerator(radius=rad, fpSize=nbs, includeChirality=True)
    return np.array(gen.GetFingerprint(mol))


def mol_to_inchikey(mol: Mol) -> str:
    """
    Converts an RDKit molecule to an InChIKey.

    :param mol: the molecule to convert
    :return: the InChIKey
    """
    return MolToInchiKey(mol)


def calc_tanimoto_similarity(arr1: NDArray[np.int8], arr2: NDArray[np.int8]) -> float:
    """
    Calculate the Tanimoto similarity between two fingerprints.

    :param arr1: the first fingerprint
    :param arr2: the second fingerprint
    :return: the Tanimoto similarity score
    """
    assert arr1.shape == arr2.shape
    assert np.all(np.logical_or(arr1 == 0, arr1 == 1))
    intersection = np.dot(arr1, arr2)
    sum_arr1 = np.sum(arr1)
    sum_arr2 = np.sum(arr2)
    score = intersection / (sum_arr1 + sum_arr2 - intersection)
    assert 0 <= score <= 1
    return score


def calc_tanimoto_similarity_rdkit(fp1: ExplicitBitVect, fp2: ExplicitBitVect) -> float:
    """
    Calculate the Tanimoto similarity between two RDKit fingerprints.

    :param fp1: the first fingerprint
    :param fp2: the second fingerprint
    :return: the Tanimoto similarity score
    """
    return TanimotoSimilarity(fp1, fp2)


def encode_mol(mol: Mol) -> str:
    """
    Encodes an RDKit molecule as a canonical isomeric SMILES string.

    :param mol: the molecule to encode
    :return: the encoded molecule
    """
    return MolToSmiles(mol, isomericSmiles=True, canonical=True)


def neutralize_mol(mol: Mol) -> None:
    """
    Neutralizes formal charges on an RDKit molecule.

    :param mol: the molecule to neutralize
    .. note:: this function modifies the input molecule in place
    """
    charge_smarts = "[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]"
    charge_pattern = smarts_to_mol(charge_smarts)
    at_matches = mol.GetSubstructMatches(charge_pattern)

    if len(at_matches) > 0:
        for match in at_matches:
            at_idx = match[0]  # get the atom index from the match tuple
            atom = mol.GetAtomWithIdx(at_idx)
            if atom.GetSymbol() in ["B"]:
                continue  # skip boron atoms
            charge = atom.GetFormalCharge()
            h_count = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(h_count - charge)
            atom.UpdatePropertyCache()


def get_default_valence(anr: int) -> int:
    """
    Returns the default valence for an atom number.

    :param anr: the atom number
    :return: the default valence
    """
    return get_periodic_table().GetDefaultValence(anr)


def get_tags_mol(mol: Mol) -> list[int]:
    """Get the atom tags from a molecule.

    :param mol: the molecule
    :return: unordered set of atom tags
    """
    tags: list[int] = []
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:
            tags.append(atom.GetIsotope())
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


def largest_fragment(mol: Mol) -> Mol:
    """
    Return the largest fragment of a molecule (by atom count).

    :param mol: input molecule
    :return: largest fragment molecule
    """
    frags: tuple[Mol, ...] = GetMolFrags(mol, asMols=True, sanitizeFrags=True)
    return max(frags, key=lambda m: m.GetNumAtoms()) if frags else mol


def count_fragments(mol: Mol) -> int:
    """
    Counts the number of fragments in a molecule.

    :param mol: the molecule to analyze
    :return: the number of fragments
    """
    frags: tuple[Mol, ...] = GetMolFrags(mol, asMols=True, sanitizeFrags=False)
    return len(frags)


def prepare_stereo(mol: Mol) -> Mol:
    """
    Ensure RDKit has the best shot at assigning stereo.
    - picks up R/S from chiral tags or 3D coordinates
    - converts wedge/hash directions to bond stereo where relevant
    - assigns E/Z on double bonds where substituent directions are known

    :param mol: input molecule
    :return: molecule with assigned stereochemistry
    """
    mol = Mol(mol)  # copy
    # If you have 2D/3D coords with wedge bonds, capture bond stereo first
    SetBondStereoFromDirections(mol)
    # Final pass to set CIP labels & E/Z where possible
    AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    return mol


def stereo_summary(mol: Mol, one_based: bool = True) -> str:
    """
    Return a compact stereochemistry summary like: "C@2=R; C@7=S; DB(3-4)=E; DB(10-11)=Z"

    :param mol: input molecule
    :param one_based: whether to use one-based atom indexing (default True)
    :return: stereochemistry summary string (none if no stereochemistry is present)
    """
    # Ensure CIP labels and bond stereo are assigned
    mol = prepare_stereo(mol)

    # Chiral centers (R/S or '?')
    chiral_bits: list[str] = []
    for idx, _ in Chem.FindMolChiralCenters(mol, includeUnassigned=True, useLegacyImplementation=False):
        # RDKit stores CIP on atom property "_CIPCode" when defined
        atom: Atom = mol.GetAtomWithIdx(idx)
        cip = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
        lbl = cip.upper() if cip else "?"
        aid = idx + 1 if one_based else idx
        chiral_bits.append(f"C@{aid}={lbl}")

    # Double bond stereo (E/Z, cis/trans, or '?')
    dbits: list[str] = []
    for b in mol.GetBonds():
        if b.GetBondType() != BondType.DOUBLE:
            continue
        st: BondStereo = b.GetStereo()
        if st in (
            BondStereo.STEREOE,
            BondStereo.STEREOZ,
            BondStereo.STEREOCIS,
            BondStereo.STEREOTRANS,
        ):
            a = b.GetBeginAtomIdx() + (1 if one_based else 0)
            c = b.GetEndAtomIdx() + (1 if one_based else 0)
            if st == BondStereo.STEREOE:
                tag = "E"
            elif st == BondStereo.STEREOZ:
                tag = "Z"
            elif st == BondStereo.STEREOCIS:
                tag = "cis"
            elif st == BondStereo.STEREOTRANS:
                tag = "trans"
            else:
                tag = "?"
            dbits.append(f"DB({a}-{c})={tag}")

    if not chiral_bits and not dbits:
        return "none"

    return ";".join(chiral_bits + dbits)


def ecfp4(mol: Mol, n_bits: int = 2048) -> ExplicitBitVect:
    """
    Compute the ECFP4 fingerprint for a molecule.

    :param mol: input molecule
    :param n_bits: number of bits in the fingerprint
    :return: ECFP4 fingerprint as an RDKit ExplicitBitVect
    """
    return GetMorganFingerprintAsBitVect(mol, radius=2, nBits=n_bits)
