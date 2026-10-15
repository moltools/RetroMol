# -*- coding: utf-8 -*-

"""Integration tests for parsing with stereochemistry matching enabled (`match_stereochemistry=True`)."""

from pathlib import Path

import pytest
from rdkit import RDLogger

from retromol.model.rules import RuleSet

from .data.integration_retromol_stereochemistry import CASES
from .helpers_retromol import assert_assemblies, load_expected_assemblies, parse_compound


# Disable RDKit warnings for cleaner test output
RDLogger.DisableLog("rdApp.*")

EXPECTED_ASSEMBLIES = load_expected_assemblies(
    Path(__file__).parent / "data" / "integration_retromol_stereochemistry_assemblies.yml"
)


def test_stereo_assembly_fixtures_match_compound_cases() -> None:
    assert set(EXPECTED_ASSEMBLIES) == {name for name, _ in CASES}


@pytest.mark.parametrize("name, smiles", CASES, ids=[c[0] for c in CASES])
def test_integration_stereochemistry_set(
    name: str,
    smiles: str,
    stereo_ruleset: RuleSet,
) -> None:
    """Graph matching preserves stereo-qualified identities such as C^Z1."""
    result = parse_compound(smiles, stereo_ruleset)
    assert_assemblies(result, EXPECTED_ASSEMBLIES[name])
