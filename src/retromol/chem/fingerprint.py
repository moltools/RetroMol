"""Module for calculating and comparing chemical fingerprints."""

from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdFingerprintGenerator import GetMorganGenerator
from rdkit.DataStructs.cDataStructs import ExplicitBitVect, TanimotoSimilarity


def mol_to_morgan_fingerprint(
    mol: Mol,
    radius: int = 2,
    num_bits: int = 2048,
    use_chirality: bool = True
) -> ExplicitBitVect:
    """
    Compute the Morgan fingerprint for a molecule.

    :param mol: input molecule
    :param radius: radius of the Morgan fingerprint
    :param num_bits: number of bits in the fingerprint
    :param use_chirality: whether to include chirality information
    :return: Morgan fingerprint as an RDKit ExplicitBitVect
    """
    generator = GetMorganGenerator(radius=radius, fpSize=num_bits, includeChirality=use_chirality)
    fingerprint = generator.GetFingerprint(mol)

    return fingerprint


def calculate_tanimoto_similarity(fp1: ExplicitBitVect, fp2: ExplicitBitVect) -> float:
    """
    Calculate the Tanimoto similarity between two RDKit fingerprints.

    :param fp1: the first fingerprint
    :param fp2: the second fingerprint
    :return: the Tanimoto similarity score
    .. note:: perfect similarity returns 1.0, no similarity returns 0.0
    """
    return TanimotoSimilarity(fp1, fp2)
