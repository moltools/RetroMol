from rdkit import Chem
import pytest

from retromol.chem import (
    Reaction,
    get_default_linearization_rules,
    get_default_sequencing_rules,
    mol_to_fingerprint,
    encode_mol,
    neutralize_mol,
)


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


def test_reaction_initialization():
    """Test initializing the Reaction class."""
    reaction_smarts = "[C:1]>>[C:1]O"
    reaction = Reaction("hydroxylation", reaction_smarts)
    assert reaction.name == "hydroxylation", "reaction name should be set correctly"
    assert reaction.reaction_smarts == reaction_smarts, "reaction SMARTS should be set correctly"
    assert reaction.rxn is not None, "reaction object should be created from SMARTS"


def test_reaction_call():
    """Test applying a Reaction to a molecule."""
    reaction_smarts = "[C:1]>>[C:1]O"
    reaction = Reaction("hydroxylation", reaction_smarts)
    reactant = Chem.MolFromSmiles("C")
    assert reactant is not None, "Reactant molecule should not be None"

    products = reaction(reactant)
    assert len(products) > 0, "reaction should produce products"
    assert len(products[0]) > 0, "Rreaction products should not be empty"


def test_reaction_with_mask():
    """Test applying a Reaction with a mask."""
    reaction_smarts = "[C:1]>>[C:1]O"
    reaction = Reaction("hydroxylation", reaction_smarts)
    reactant = Chem.MolFromSmiles("[13C]C")
    assert reactant is not None, "reactant molecule should not be None"

    # apply reactionwith a mask; without setting unique tags, the reaction will throw ValueError
    with pytest.raises(ValueError):
        products = reaction(reactant, mask={13})

    # apply reaction with a mask; make sure all atoms have unique tags
    taken_tags = [atom.GetIsotope() for atom in reactant.GetAtoms() if atom.GetIsotope() != 0]
    for atom in reactant.GetAtoms():
        if atom.GetIsotope() == 0:
            tag = 1
            while tag in taken_tags:
                tag += 1
            atom.SetIsotope(tag)
            taken_tags.append(tag)
    products = reaction(reactant, mask={13})
    assert len(products) > 0, "reaction with mask should produce products"


def test_get_default_linearization_rules():
    """Test getting default linearization rules."""
    rules = get_default_linearization_rules()
    assert isinstance(rules, list), "linearization rules should be a list"
    assert all(isinstance(rule, Reaction) for rule in rules), "linearization rules should contain Reaction objects"


def test_get_default_sequencing_rules():
    """Test getting default sequencing rules."""
    rules = get_default_sequencing_rules()
    assert isinstance(rules, list), "sequencing rules should be a list"
    assert all(isinstance(rule, Reaction) for rule in rules), "sequencing rules should contain Reaction objects"
