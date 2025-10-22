# -*- coding: utf-8 -*-

"""Pytest configuration for loading rule sets and wave configurations."""

import pytest

from .helpers import load_rule_set, load_wave_config


@pytest.fixture(scope="session")
def rule_set():
    """Load rule set once per test session."""
    return load_rule_set()


@pytest.fixture(scope="session")
def wave_config():
    """Load a simple wave configuration once per test session."""
    return load_wave_config()
