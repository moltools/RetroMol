# -*- coding: utf-8 -*-

"""Integration tests for the demo set of compounds."""

import pytest
from rdkit import RDLogger

from retromol.model.rules import RuleSet

from .data.integration_demo_set import CASES
from .helpers import assert_result, parse_compound


# Disable RDKit warnings for cleaner test output
RDLogger.DisableLog("rdApp.*")


@pytest.mark.parametrize("name, smiles, expected_coverage, expected_monomers", CASES, ids=[c[0] for c in CASES])
def test_integration_demo_set(name: str, smiles: str, expected_coverage: float, expected_monomers: list[list[str]], ruleset: RuleSet) -> None:
    """
    Integration test for the demo set of compounds.

    :param name: the name of the test case
    :param smiles: the SMILES string of the compound to test
    :param expected_coverage: the expected total coverage value
    :param expected_monomers: the expected list of monomer identities
    :param ruleset: the RuleSet to use for parsing
    """
    print(f"testing {name}...")
    result = parse_compound(smiles, ruleset)
    assert_result(result, expected_coverage, expected_monomers)
