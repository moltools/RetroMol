from rdkit import Chem

from retromol.chem import mol_to_fingerprint, encode_mol, neutralize_mol


def test_mol_to_fingerprint():
    """Test the mol_to_fingerprint function."""
    mol = Chem.MolFromSmiles("CCO")  # Ethanol molecule
    assert mol is not None, "molecule should not be None"

    fp = mol_to_fingerprint(mol, radius=2, num_bits=2048)
    assert fp.GetNumBits() == 2048, "fingerprint should have the correct number of bits"


def test_encode_mol():
    """Test the encode_mol function."""
    mol = Chem.MolFromSmiles("CCO")  # Ethanol molecule
    assert mol is not None, "molecule should not be None"

    encoded = encode_mol(mol)
    assert isinstance(encoded, int), "encoded molecule should be an integer"

    # test with additional_bits
    additional_bits = [1, 0, 1]
    encoded_with_bits = encode_mol(mol, additional_bits=additional_bits)
    assert isinstance(encoded_with_bits, int), "encoded molecule with additional bits should be an integer"
    assert encoded_with_bits != encoded, "encoding with additional bits should change the result"


def test_neutralize_mol():
    """Test the neutralize_mol function."""
    # Molecule with a charged group (ammonium ion)
    mol = Chem.MolFromSmiles("[NH4+]")
    assert mol is not None, "molecule should not be None"

    # neutralize the molecule
    neutralize_mol(mol)
    atom = mol.GetAtomWithIdx(0)
    assert atom.GetFormalCharge() == 0, "atom should be neutralized"
    assert atom.GetTotalNumHs() == 3, "atom should have the correct number of hydrogens"

    # test with a neutral molecule (should remain unchanged)
    neutral_mol = Chem.MolFromSmiles("CCO")
    neutralize_mol(neutral_mol)
    assert neutral_mol.GetNumAtoms() == 3, "neutral molecule should remain unchanged"
