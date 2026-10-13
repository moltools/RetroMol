# -*- coding: utf-8 -*-

"""Compare parsed assemblies with explicit graph expectations for demo compounds."""

from pathlib import Path

import pytest
from rdkit import RDLogger

from retromol.model.rules import RuleSet

from .data.integration_retromol import CASES
from .helpers_retromol import assert_assemblies, load_expected_assemblies, parse_compound


# Disable RDKit warnings for cleaner test output
RDLogger.DisableLog("rdApp.*")

EXPECTED_ASSEMBLIES = load_expected_assemblies(
    Path(__file__).parent / "data" / "integration_retromol_assemblies.yml"
)


def test_assembly_fixtures_match_compound_cases() -> None:
    assert set(EXPECTED_ASSEMBLIES) == {name for name, _ in CASES}


@pytest.mark.parametrize("name, smiles", CASES, ids=[c[0] for c in CASES])
def test_integration_demo_set(name: str, smiles: str, ruleset: RuleSet) -> None:
    """Every expected graph must occur among the valid parsed assemblies."""
    result = parse_compound(smiles, ruleset)
    assert_assemblies(result, EXPECTED_ASSEMBLIES[name])
