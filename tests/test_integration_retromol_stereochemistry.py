# -*- coding: utf-8 -*-

"""Integration tests for parsing with stereochemistry matching enabled (`match_stereochemistry=True`)."""

import pytest
from rdkit import RDLogger

from retromol.model.rules import RuleSet

from .data.integration_retromol_stereochemistry import CASES
from .helpers_retromol import assert_result, parse_compound


# Disable RDKit warnings for cleaner test output
RDLogger.DisableLog("rdApp.*")


@pytest.mark.parametrize("name, smiles, expected_coverage, expected_monomers", CASES, ids=[c[0] for c in CASES])
def test_integration_stereochemistry_set(
    name: str,
    smiles: str,
    expected_coverage: float,
    expected_monomers: list[str],
    stereo_ruleset: RuleSet,
) -> None:
    """
    Integration test for compounds parsed with stereochemistry matching enabled.

    :param name: the name of the test case
    :param smiles: the SMILES string of the compound to test
    :param expected_coverage: the expected total coverage value
    :param expected_monomers: the expected list of monomer identities
    :param stereo_ruleset: the RuleSet (with match_stereochemistry=True) to use for parsing
    """
    print(f"testing {name}...")
    result = parse_compound(smiles, stereo_ruleset)
    assert_result(result, expected_coverage, expected_monomers)
