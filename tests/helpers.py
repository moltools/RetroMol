"""Shared helpers for RetroMol integration tests."""

from retromol.model.rules import RuleSet
from retromol.model.result import Result
from retromol.model.submission import Submission
from retromol.model.reaction_graph import MolNode
from retromol.pipelines.parsing import run_retromol_with_timeout


def load_rule_set() -> RuleSet:
    """
    Load the default RetroMol rule set once.
    
    :return: the loaded RuleSet object
    """
    return RuleSet.load_default(match_stereochemistry=False)


def parse_compound(smiles: str, ruleset: RuleSet) -> Result:
    """
    Parse a compound SMILES string into an Result object.

    :param smiles: the SMILES string of the compound to parse
    :param ruleset: the RuleSet to use for parsing
    :return: the resulting Result object
    """
    submission = Submission(smiles)
    return run_retromol_with_timeout(submission, ruleset)


def compare_floats(a: float, b: float, tol: float = 1e-2) -> bool:
    """
    Compare two floating-point numbers for equality within a tolerance.

    :param a: the first float to compare
    :param b: the second float to compare
    :param tol: the tolerance for comparison
    :return: True if the numbers are equal within the tolerance, False otherwise
    """
    return abs(a - b) <= tol


def compare_lists(list1: list[str], list2: list[str]) -> bool:
    """
    Compare two lists of strings for equality, ignoring order.

    :param list1: the first list to compare
    :param list2: the second list to compare
    :return: True if the lists contain the same elements, False otherwise
    """
    return sorted(list1) == sorted(list2)


def assert_result(result: Result, expected_coverage: float, expected_monomers: list[str]) -> None:
    """
    Common assertion logic used by all integration tests.
    
    :param result: the Result object to check
    :param expected_coverage: the expected total coverage value
    :param expected_monomers: the expected list of monomer identities
    """
    coverage: float = result.calculate_coverage()
    assert compare_floats(coverage, expected_coverage), f"expected coverage {expected_coverage}, got {coverage}"

    ident_nodes: MolNode = result.reaction_graph.identified_nodes.values()
    assert all(n.is_identified for n in ident_nodes), "not all identified nodes are marked as identified"
    found_monomers: list[str] = [n.identity.name for n in ident_nodes]

    # Sort monomers before comparison; easier to read in case of failure
    found_monomers.sort()
    expected_monomers.sort()

    assert compare_lists(found_monomers, expected_monomers), f"expected monomers {expected_monomers}, got {found_monomers}"
