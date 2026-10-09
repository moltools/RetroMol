"""Pytest configuration for loading rule sets and wave configurations."""

import pytest

from retromol.model.rules import RuleSet

from .helpers_retromol import load_rule_set, load_stereo_rule_set


@pytest.fixture(scope="session")
def ruleset() -> RuleSet:
    """
    Load rule set once per test session.

    :return: the loaded RuleSet object
    """
    return load_rule_set()


@pytest.fixture(scope="session")
def stereo_ruleset() -> RuleSet:
    """
    Load rule set once per test session, matching stereochemistry.

    :return: the loaded RuleSet object
    """
    return load_stereo_rule_set()
