# -*- coding: utf-8 -*-

"""Compare demo compounds with fixed graph and ordered readout expectations."""

from pathlib import Path

import pytest
from rdkit import RDLogger

from retromol.model.rules import RuleSet

from .data.integration_retromol import CASES
from .helpers_retromol import (
    assert_assemblies, assert_linear_readout, load_expected_assemblies,
    load_expected_readouts, parse_compound,
)


# Disable RDKit warnings for cleaner test output
RDLogger.DisableLog("rdApp.*")

EXPECTED_ASSEMBLIES = load_expected_assemblies(
    Path(__file__).parent / "data" / "integration_retromol_assemblies.yml"
)
EXPECTED_READOUTS = load_expected_readouts(
    Path(__file__).parent / "data" / "integration_retromol_readouts.yml", EXPECTED_ASSEMBLIES,
)


def test_assembly_fixtures_match_compound_cases() -> None:
    assert set(EXPECTED_ASSEMBLIES) == {name for name, _ in CASES}
    assert set(EXPECTED_READOUTS) == set(EXPECTED_ASSEMBLIES)


@pytest.mark.parametrize("name, smiles", CASES, ids=[c[0] for c in CASES])
def test_integration_demo_set(name: str, smiles: str, ruleset: RuleSet) -> None:
    """Every required graph must have its sequence; also check the default view."""
    result = parse_compound(smiles, ruleset)
    readouts = EXPECTED_READOUTS[name]
    assert_assemblies(result, EXPECTED_ASSEMBLIES[name], expected_readouts=readouts.assemblies)
    assert_linear_readout(result.assembly_graph, readouts.default)
    assert result.linear_readout.primary_sequence_parts() == [list(part) for part in readouts.default]


@pytest.mark.parametrize("name, smiles", [
    pytest.param(name, smiles, id=name, marks=pytest.mark.xfail(
        strict=True, raises=AssertionError, reason=EXPECTED_READOUTS[name].known_starter_order_issue,
    ))
    for name, smiles in CASES if EXPECTED_READOUTS[name].known_starter_order_issue
])
def test_known_backbone_starter_order(name, smiles, ruleset):
    """Keep known ordering errors visible independently of compatibility snapshots.

    These checks express the desired property, not the current faulty order.
    A fix must remove its expected-failure mark and update the corresponding
    regression fixture; strict XPASS prevents leaving a fixed issue hidden.
    """
    result = parse_compound(smiles, ruleset)
    paths = result.linear_readout.primary_sequence_parts()
    single_starter_paths = [part for part in paths if sum(token.endswith(' starter') for token in part) == 1]
    assert single_starter_paths, f"{name}: the starter disappeared from the readout"
    assert all(part[0].endswith(' starter') for part in single_starter_paths), paths
