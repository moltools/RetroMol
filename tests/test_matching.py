from rdkit import Chem
import pytest

from retromol.matching import Motif, get_default_motifs


def test_motif_initialization():
    """Test initializing the Motif class."""
    motif = Motif("hydroxyl", "O")
    assert motif.name == "hydroxyl", "motif name should be set correctly"
    assert motif.smiles == "O", "motif SMILES should be set correctly"
    assert motif.mol is not None, "motif molecule should be created from SMILES"


def test_motif_is_match():
    """Test the is_match method of Motif."""
    motif = Motif("hydroxyl", "O")
    mol = Chem.MolFromSmiles("O")
    assert motif.is_match(mol), "motif should match the molecule"

    non_matching_mol = Chem.MolFromSmiles("CC")
    assert not motif.is_match(non_matching_mol), "motif should not match a non-matching molecule"


def test_get_default_motifs():
    """Test getting default motifs."""
    motifs = get_default_motifs()
    assert isinstance(motifs, list), "default motifs should be a list"
    assert all(isinstance(motif, Motif) for motif in motifs), "default motifs should contain Motif objects"
