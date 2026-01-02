"""Module defining reaction and matching rules."""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from importlib.resources import files
from pathlib import Path
from typing import Any

import yaml
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdChemReactions import ChemicalReaction

import retromol.data
from retromol.chem.mol import (
    mol_to_smiles,
    smiles_to_mol,
    count_fragments,
    sanitize_mol,
    reassign_stereochemistry,
)
from retromol.chem.reaction import smarts_to_reaction
from retromol.chem.tagging import get_tags_mol, remove_tags
from retromol.chem.masking import is_masked_preserved


log = logging.getLogger(__name__)


@dataclass(frozen=True)
class ReactionRule:
    """
    Represents a chemical reaction rule defined by a SMARTS pattern.

    :var name: str: name of the reaction rule
    :var smarts: str: SMARTS pattern defining the reaction
    :var props: dict[str, Any]: additional properties associated with the rule
    """

    name: str
    smarts: str
    props: dict[str, Any]

    rxn: ChemicalReaction = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """
        Initialize the ChemicalReaction from the SMARTS pattern.
        """
        rxn = smarts_to_reaction(self.smarts)
        object.__setattr__(self, "rxn", rxn)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ReactionRule":
        """
        Create a ReactionRule instance from a dictionary.

        :param data: dict[str, Any]: dictionary containing rule data
        :return: ReactionRule: the created ReactionRule instance
        """
        return ReactionRule(
            name=data["name"],
            smarts=data["smarts"],
            props=data.get("props", {}),
        )
    
    def apply(self, reactant: Mol, mask_tags: set[int] | None = None) -> list[list[Mol]]:
        """
        Apply the reaction to the given reactant molecule, optionally enforcing a mask on atom tags.

        :param reactant: Mol: the reactant molecule
        :param mask_tags: set[int] | None: set of atom tags (isotope-based tags) that are allowed to change
        :return: list[list[Mol]]: list of unique product tuples (each tuple as a list[Mol])
        """ 
        log.debug(f"applying reaction rule '{self.name}'")

        results = self.rxn.RunReactants([reactant])
        if not results:
            log.debug("no valid products generated for reactant")
            return []
        
        # Sanitize and filter
        kept: list[list[Mol]] = []
        for tup in results:
            products: list[Mol] = []

            # Quick shape check, and sanitize
            atom_tag_sets: list[set[int]] = []
            ok_tuple = True
            for prod in tup:

                # Check if product is single component
                if not count_fragments(prod) == 1:
                    log.debug("product has multiple components, skipping")
                    ok_tuple = False
                    break
                
                # Sanitize in place
                if not sanitize_mol(prod, correct_hydrogens=True):
                    log.debug("product sanitization failed, skipping")
                    ok_tuple = False
                    break

            # Reassign stereo on the sanitized product
            prod = reassign_stereochemistry(prod)

            products.append(prod)
            atom_tag_sets.append(get_tags_mol(prod))

            if not ok_tuple:
                log.debug("product tuple failed validation, skipping")
                continue

            # Disallow overlapping tag sets across products
            total_tags = sum(len(s) for s in atom_tag_sets)
            union_tags = len(set().union(*atom_tag_sets)) if atom_tag_sets else 0
            if atom_tag_sets and total_tags != union_tags:
                log.debug("products share atom tags, skipping")
                continue

            # Mask check
            if mask_tags is not None and not is_masked_preserved(reactant, products, mask_tags):
                log.debug("products modify tags outside mask, skipping")
                continue

            kept.append(products)

        if len(kept) <= 1:
            return kept
        
        # Stereo-aware dereplication (order-insensitive, multiplicity-preserving)
        seen: dict[tuple[tuple[str, int], ...], int] = {}
        unique: list[list[Mol]] = []
        for res in kept:

            # Create keys based on the SMILES of products without tags
            c = Counter(mol_to_smiles(p, include_tags=False, isomeric=True, canonical=True) for p in res)
            key = tuple(sorted(c.items(), key=lambda kv: kv[0]))

            if key in seen:
                continue

            seen[key] = 1
            unique.append(res)
        
        return unique


@dataclass(frozen=True)
class MatchingRule:
    """
    Represents a molecular matching rule defined by a SMILES pattern.

    :var name: str: name of the matching rule
    :var smiles: str: SMILES pattern defining the motif
    :var props: dict[str, Any]: additional properties associated with the rule
    """

    name: str
    smiles: str
    props: dict[str, Any]

    mol: Mol = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """
        Initialize the Mol from the SMILES pattern.
        """
        mol = smiles_to_mol(self.smiles)
        object.__setattr__(self, "mol", mol)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "MatchingRule":
        """
        Create a MatchingRule instance from a dictionary.

        :param data: dict[str, Any]: dictionary containing rule data
        :return: MatchingRule: the created MatchingRule instance
        """
        return MatchingRule(
            name=data["name"],
            smiles=data["smiles"],
            props=data.get("props", {}),
        )

    def is_match(self, mol: Mol, match_stereochemistry: bool = False) -> bool:
        """
        Check if the given molecule matches this rule.

        :param mol: Mol: molecule to check
        :param match_stereochemistry: bool: whether to consider stereochemistry in matching
        :return: bool: True if the molecule matches the rule, False otherwise
        """
        has_substruct_match = mol.HasSubstructMatch(self.mol, useChirality=match_stereochemistry)
        has_equal_num_atoms = mol.GetNumAtoms() == self.mol.GetNumAtoms()
        has_equal_num_bonds = mol.GetNumBonds() == self.mol.GetNumBonds()
        
        if has_substruct_match and has_equal_num_atoms and has_equal_num_bonds:
            return True

        return False


@dataclass(frozen=True)
class RuleSet:
    """
    Represents a set of reaction and matching rules.

    :var reaction_rules: list[ReactionRule]: list of reaction rules
    :var matching_rules: list[MatchingRule]: list of matching rules
    """

    reaction_rules: list[ReactionRule]
    matching_rules: list[MatchingRule]

    @classmethod
    def load_default(cls) -> "RuleSet":
        """
        Load the default set of reaction and matching rules.

        :return: RuleSet: the default rule set
        """
        path_reaction_rules = Path(files(retromol.data).joinpath("default_reaction_rules.yml"))
        path_matching_rules = Path(files(retromol.data).joinpath("default_matching_rules.yml"))

        with open(path_reaction_rules, "r") as fo:
            reaction_rules_data = yaml.safe_load(fo)

        with open(path_matching_rules, "r") as fo:
            matching_rules_data = yaml.safe_load(fo)

        reaction_rules = [ReactionRule.from_dict(d) for d in reaction_rules_data]
        matching_rules = [MatchingRule.from_dict(d) for d in matching_rules_data]

        return RuleSet(reaction_rules, matching_rules)
