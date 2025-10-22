# -*- coding: utf-8 -*-

"""Shared helpers for RetroMol integration tests."""

from __future__ import annotations

from importlib.resources import files
from typing import Any, Dict, List

import yaml

import retromol.data
from retromol import api, io, readout, rules


def load_rule_set() -> rules.Rules:
    """Load the default RetroMol rule set once."""
    path_reaction_rules = str(files(retromol.data).joinpath("default_reaction_rules.yml"))
    path_matching_rules = str(files(retromol.data).joinpath("default_matching_rules.yml"))
    return rules.load_rules_from_files(path_reaction_rules, path_matching_rules)


def load_wave_config() -> Dict[str, Any]:
    """Load the default wave configuration once."""
    path_wave_config = str(files(retromol.data).joinpath("default_wave_config.yml"))

    with open(path_wave_config) as f:
        return yaml.safe_load(f)


def parse_compound(
    smiles: str,
    rule_set: rules.Rules,
    wave_config: Dict[str, Any],
    *,
    match_stereochemistry: bool = False,
) -> io.Result:
    """Parse a compound SMILES string into an io.Result object."""
    mol = io.Input("test_compound", smiles)
    return api.run_retromol_with_timeout(mol, rule_set, wave_config, match_stereochemistry=match_stereochemistry)


def compare_lists_of_lists(a: List[List[str]], b: List[List[str]]) -> bool:
    # Convert each inner list to a frozenset (hashable, unordered)
    set_a = {frozenset(inner) for inner in a}
    set_b = {frozenset(inner) for inner in b}
    return set_a == set_b


def assert_result(result: io.Result, expected_coverage: float, expected_mappings: List[List[str]]) -> None:
    """Common assertion logic used by all integration tests."""
    best_total_coverage: float = result.best_total_coverage()
    assert best_total_coverage == expected_coverage, f"Expected coverage {expected_coverage}, got {best_total_coverage}"

    mappings = readout.optimal_mappings_with_timeout(result)
    parsed_mappings: List[List[str]] = []
    for mapping in mappings:
        parsed_mapping: List[str] = []
        for node in mapping["nodes"]:
            parsed_mapping.append(node["identity"])
        parsed_mapping.sort()
        parsed_mappings.append(parsed_mapping)

    assert compare_lists_of_lists(expected_mappings, parsed_mappings), (
        f"Expected mappings {expected_mappings}, got {parsed_mappings}"
    )
