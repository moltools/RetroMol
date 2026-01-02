"""Data structure for representing a molecular identity."""

from dataclasses import dataclass

from retromol.model.rules import MatchingRule


@dataclass(frozen=True)
class MolIdentity:
    """
    Represents the identity of a molecule based on matched rules.

    :var matched_rules: list[str]: List of matched rule identifiers
    """

    matched_rule: MatchingRule

    @property
    def name(self) -> str:
        """
        Get the name of the matched rule.

        :return: name of the matched rule
        """
        return self.matched_rule.name
