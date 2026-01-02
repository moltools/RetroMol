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

    :param mol: RDKit molecule to match
    :param rules: list of MatchingRule to use for matching
    :param match_stereochemistry: whether to consider stereochemistry in matching
    :return: MolIdentity | None: the identity if matched, else None
    .. note:: this function uses a greedy approach to match a molecule to a motif
    """
    for rl in rules:
        if rl.is_match(mol, match_stereochemistry):
            return MolIdentity(matched_rule=rl)

    return None
