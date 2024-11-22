from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


def mol_to_fingerprint(mol: Chem.Mol, radius: int = 2, num_bits: int = 2048) -> rdFingerprintGenerator.ExplicitBitVect:
    """Converts an RDKit molecule to a Morgan fingerprint.

    :param mol: The molecule to convert.
    :type mol: Chem.Mol
    :param radius: The radius of the fingerprint.
    :type radius: int
    :param num_bits: The number of bits in the fingerprint.
    :type num_bits: int
    :return: The Morgan fingerprint.
    :rtype: ExplicitBitVect
    """
    fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=num_bits)
    return fp_generator.GetFingerprint(mol)


def encode_mol(mol: Chem.Mol) -> int:
    """Encodes an RDKit molecule as an integer.

    :param mol: The molecule to encode.
    :type mol: Chem.Mol
    :return: The encoded molecule.
    :rtype: int
    """
    fp = mol_to_fingerprint(mol)
    fp_hash = hash(fp.data.tobytes())
    return fp_hash


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
