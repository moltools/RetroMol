# -*- coding: utf-8 -*-

"""Integration test cases for parsing with stereochemistry matching enabled (`match_stereochemistry=True`).

Unlike `integration_retromol.py`, these cases assert on stereo-qualified monomer
identities (e.g. `C^Z1` instead of `C1`) and exist to catch regressions in how E/Z and
R/S stereo descriptors are carried through the reaction graph, in particular double
bond stereo that is only exposed once a ring is opened mid-pipeline (see the
"stereo registry" mechanism in `retromol.chem.stereo`).
"""

# Cases are formatted as:
#
# (
#    name,
#    smiles,
#    expected coverage_score,
#    list of found monomers
# ),

CASES: list[tuple[str, str, float, list[str]]] = [
    (
        "ring-opened Z double bond (small-ring C1 unit)",
        r"O[C@@H](CCC)C[C@H](C[C@@H](C[C@H](C[C@H](CC=C1)OC1=O)O)C)OC",
        1.0,
        ["butanoic acid", "acetic acid", "D1", "B^S1", "B^R1", "D^R8", "B^R1", "B^R1", "C^Z1", "methylation"],
    ),
]
