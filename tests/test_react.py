from rdkit import Chem
import pytest

from retromol.chem import (
    Motif,
    Reaction,
    get_default_motifs,
    get_default_linearization_rules,
    get_default_sequencing_rules,
)

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


def test_get_default_motifs():
    """Test getting default motifs."""
    motifs = get_default_motifs()
    assert isinstance(motifs, list), "default motifs should be a list"
    assert all(isinstance(motif, Motif) for motif in motifs), "default motifs should contain Motif objects"


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
