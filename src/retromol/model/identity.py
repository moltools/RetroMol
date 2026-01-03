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
    
    @property
    def terminal(self) -> bool:
        """
        Check if the matched rule indicates a terminal identity.

        :return: True if the matched rule is terminal, False otherwise
        """
        return self.matched_rule.terminal
    
    def to_dict(self) -> dict:
        """
        Serialize the MolIdentity to a dictionary.

        :return: dictionary representation of the MolIdentity
        """
        return {
            "matched_rule": self.matched_rule.to_dict(),
        }
    
    @classmethod
    def from_dict(cls, data: dict) -> "MolIdentity":
        """
        Deserialize a MolIdentity from a dictionary.

        :param data: dictionary representation of the MolIdentity
        :return: MolIdentity object
        """
        matched_rule = MatchingRule.from_dict(data["matched_rule"])
        return cls(
            matched_rule=matched_rule,
        )
