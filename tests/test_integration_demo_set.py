# -*- coding: utf-8 -*-

"""Integration tests for the demo set of compounds."""

from typing import Any, Dict, List

import pytest

from retromol import rules

from .data.integration_demo_set import CASES
from .helpers import assert_result, parse_compound


@pytest.mark.parametrize("identifier, smiles, expected_coverage, expected_mappings", CASES, ids=[c[0] for c in CASES])
def test_integration_demo_set(
    identifier: str,
    smiles: str,
    expected_coverage: float,
    expected_mappings: List[List[str]],
    rule_set: rules.Rules,
    wave_config: Dict[str, Any],
) -> None:
    print(f"Testing {identifier}...")
    result = parse_compound(smiles, rule_set, wave_config, match_stereochemistry=False)
    assert_result(result, expected_coverage, expected_mappings)
