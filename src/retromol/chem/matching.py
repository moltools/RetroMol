"""Module for matching molecules to MatchingRules."""

from retromol.chem.mol import Mol
from retromol.model.rules import MatchingRule
from retromol.model.identity import MolIdentity


def match_mol(
    mol: Mol,
    rules: list[MatchingRule],
    match_stereochemistry: bool = False
) -> MolIdentity | None:
    """
    Match a molecule to a motif.

    :param mol: RDKit molecule to match.
    :param rules: List of MatchingRule to use for matching.
    :param match_stereochemistry: Whether to consider stereochemistry in matching.
    :return: MolIdentity | None: The identity if matched, else None.
    .. note:: This function uses a greedy approach to match a molecule to a motif.
    """
    for rl in rules:
        if not rl.stereochemistry or (rl.stereochemistry and match_stereochemistry):
            resolved = rl.match(mol, match_stereochemistry)
            if resolved is not None:
                return MolIdentity(matched_rule=resolved)

    return None
