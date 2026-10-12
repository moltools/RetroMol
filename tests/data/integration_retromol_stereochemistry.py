# -*- coding: utf-8 -*-

"""Compound inputs for stereo-aware assembly-graph integration tests.

Required assembly graphs, including coverage, live in
integration_retromol_stereochemistry_assemblies.yml. Cases contain (name, SMILES).
"""

CASES: list[tuple[str, str]] = [
    (
        "ring-opened Z double bond (small-ring C1 unit)",
        r"O[C@@H](CCC)C[C@H](C[C@@H](C[C@H](C[C@H](CC=C1)OC1=O)O)C)OC",
    ),
]
