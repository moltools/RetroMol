# -*- coding: utf-8 -*-

"""Integration tests for parsing with stereochemistry matching enabled (`match_stereochemistry=True`)."""

from pathlib import Path

import pytest
from rdkit import RDLogger

from retromol.model.rules import RuleSet

from .data.integration_retromol_stereochemistry import CASES
from .helpers_retromol import (
    assert_assemblies, assert_linear_readout, load_expected_assemblies,
    load_expected_readouts, parse_compound,
)


# Disable RDKit warnings for cleaner test output
RDLogger.DisableLog("rdApp.*")

EXPECTED_ASSEMBLIES = load_expected_assemblies(
    Path(__file__).parent / "data" / "integration_retromol_stereochemistry_assemblies.yml"
)
EXPECTED_READOUTS = load_expected_readouts(
    Path(__file__).parent / "data" / "integration_retromol_stereochemistry_readouts.yml", EXPECTED_ASSEMBLIES,
)


def test_stereo_assembly_fixtures_match_compound_cases() -> None:
    assert set(EXPECTED_ASSEMBLIES) == {name for name, _ in CASES}
    assert set(EXPECTED_READOUTS) == set(EXPECTED_ASSEMBLIES)


@pytest.mark.parametrize("name, smiles", CASES, ids=[c[0] for c in CASES])
def test_integration_stereochemistry_set(
    name: str,
    smiles: str,
    stereo_ruleset: RuleSet,
) -> None:
    """Graphs and ordered sequences preserve stereo-qualified identities such as C^Z1."""
    result = parse_compound(smiles, stereo_ruleset)
    readouts = EXPECTED_READOUTS[name]
    assert_assemblies(result, EXPECTED_ASSEMBLIES[name], expected_readouts=readouts.assemblies)
    assert_linear_readout(result.assembly_graph, readouts.default)
    assert result.linear_readout.primary_sequence_parts() == [list(part) for part in readouts.default]
