"""Module defining reaction and matching rules."""

import logging
import itertools
import hashlib
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
from retromol.chem.reaction import smarts_to_reaction, reactive_template_atoms
from retromol.chem.tagging import get_tags_mol
from retromol.chem.masking import is_masked_preserved


log = logging.getLogger(__name__)


@dataclass(frozen=True)
class ReactionRule:
    """
    Represents a chemical reaction rule defined by a SMARTS pattern.

    :var name: str: name of the reaction rule
    :var smarts: str: SMARTS pattern defining the reaction
    :var props: dict[str, Any]: additional properties associated with the rule
    :var allowed_in_bulk: bool: whether this rule is allowed to be applied in bulk preprocessing
    """

    name: str
    smarts: str
    props: dict[str, Any]
    allowed_in_bulk: bool = False

    rxn: ChemicalReaction = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """
        Initialize the ChemicalReaction from the SMARTS pattern.
        """
        rxn = smarts_to_reaction(self.smarts)
        object.__setattr__(self, "rxn", rxn)

    @property
    def id(self) -> str:
        """
        Unique identifier for the reaction rule based on its SMARTS pattern.

        :return: str: unique identifier
        """
        return hashlib.sha256(self.smarts.encode("utf-8")).hexdigest()

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ReactionRule":
        """
        Create a ReactionRule instance from a dictionary.

        :param data: dict[str, Any]: dictionary containing rule data
        :return: ReactionRule: the created ReactionRule instance
        """
        reaction_rule = cls(
            name=data["name"],
            smarts=data["smarts"],
            props=data.get("props", {}),
            allowed_in_bulk=data.get("allowed_in_bulk", False),
        )
        return reaction_rule
    
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
            log.debug("no products generated for reactant")
            return []
        log.debug(f"generated {len(results)} raw product tuple(s)")
        
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
                if not sanitize_mol(prod, fix_hydrogens=True):
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


def index_uncontested(
    mol: Mol,
    rules: list[ReactionRule],
    failed_combos: set[tuple[int, frozenset[int]]],
) -> list[tuple[ReactionRule, set[int]]]:
    """
    Index uncontested reactions for applying preprocessing rules in bulk.

    :param mol: RDKit molecule
    :param rules: List of preprocessing rules
    :param failed_combos: Set of failed combinations to avoid infinite loops
    :return: Uncontested reactions
    """
    up_for_election: list[tuple[ReactionRule, set[int], set[int]]] = []
    for rl in rules:
        if not rl.rxn:
            continue  # skip rules without a reaction template

        reactive_inds = reactive_template_atoms(rl.rxn)[0]
        all_reactant_matches: list[tuple[tuple[int, ...], ...]] = []
        all_reactant_matches_reactive_items: list[list[list[int]]] = []
        for template_idx in range(rl.rxn.GetNumReactantTemplates()):
            reactant_template = rl.rxn.GetReactantTemplate(template_idx)
            reactant_matches: tuple[tuple[int, ...], ...] = mol.GetSubstructMatches(reactant_template)
            all_reactant_matches.append(reactant_matches)
            new_reactant_matches: list[list[int]] = []
            for reactant_match in reactant_matches:
                new_reactant_matches.append([reactant_match[idx] for idx in reactive_inds])
            all_reactant_matches_reactive_items.append(new_reactant_matches)

        # Generate all possible match sets, for when reaction template matches multiple sites
        match_sets = list(itertools.product(*all_reactant_matches))
        match_sets_reactive_items = list(itertools.product(*all_reactant_matches_reactive_items))
        match_sets = [set(itertools.chain(*match_set)) for match_set in match_sets]
        match_sets_reactive_items = [set(itertools.chain(*match_set)) for match_set in match_sets_reactive_items]
        for match_set, match_set_reactive_items in zip(match_sets, match_sets_reactive_items, strict=True):
            up_for_election.append((rl, match_set, match_set_reactive_items))

    # Check which reactions with matched templates are uncontested and which are contested
    uncontested: list[tuple[ReactionRule, set[int]]] = []
    for i, (rl, match_set, match_set_reactive_items) in enumerate(up_for_election):
        # TODO: Rules with ring matching conditions are always contested
        # if rl.has_ring_matching_condition():
        #     continue

        # Check if match set has overlap with any other match set
        # has_overlap = any(match_set.intersection(o) for j, (_, o, o_r) in enumerate(up_for_election) if i != j)
        has_overlap = any(match_set_reactive_items.intersection(o_r) for j, (_, _, o_r) in enumerate(up_for_election) if i != j)
        if not has_overlap:
            uncontested.append((rl, match_set))

    # Filter out failed combinations to avoid infinite loops
    uncontested = [(rl, match_set) for rl, match_set in uncontested if (rl.id, frozenset(match_set)) not in failed_combos]

    return uncontested


def apply_uncontested(
    parent: Mol,
    uncontested: list[tuple[ReactionRule, set[int]]],
    original_taken_tags: set[int],
) -> tuple[list[Mol], list[tuple[ReactionRule, set[int]]], set[tuple[int, frozenset[int]]]]:
    """
    Apply uncontested reactions in bulk.

    :param parent: RDKit molecule
    :param uncontested: List of uncontested reactions
    :param original_taken_tags: List of atom tags from original reactant
    :return: list of true products, a list of applied ReactionRules with their masks,  and a set of failed combinations
    """
    applied_reactions: list[tuple[ReactionRule, set[int]]] = []

    tags_in_parent: set[int] = set(get_tags_mol(parent))

    # We make sure all atoms, even the ones not from original reactant, have a
    # unique isotope number, so we can track them through consecutive reactions
    temp_taken_tags = get_tags_mol(parent)
    for atom in parent.GetAtoms():
        if atom.GetIsotope() == 0:
            tag = 1
            while tag in original_taken_tags or tag in temp_taken_tags:
                tag += 1
            atom.SetIsotope(tag)
            temp_taken_tags.add(tag)

    # Validate that all atoms have a unique tag
    num_tagged_atoms = len(set(get_tags_mol(parent)))
    if num_tagged_atoms != len(parent.GetAtoms()):
        raise ValueError("Not all atoms have a unique tag before applying uncontested reactions")

    # Map tags to atomic nums so we can create masks and reassign atomic nums later on
    idx_to_tag = {a.GetIdx(): a.GetIsotope() for a in parent.GetAtoms()}

    # All uncontested reactions become a single node in the reaction_graph
    products: list[Mol] = []
    failed_combos: set[tuple[int, frozenset[int]]] = set()  # keep track of failed combinations to avoid infinite loops

    for rl, match_set in uncontested:
        msk = set([idx_to_tag[idx] for idx in match_set])  # create mask for reaction

        # We use the input parent if there are no products, otherwise we have to find out
        # which product now contains the mask (i.e., the reaction template for this reaction)
        if len(products) != 0:
            new_parent: Mol | None = None
            for product in products:
                product_tags = set(get_tags_mol(product))
                if msk.issubset(product_tags):
                    new_parent = product
                    products = [p for p in products if p != product]
                    break

            if new_parent is None:
                # raise ValueError("no product found that contains the mask")
                # If no product is found, we continue with the next uncontested reaction
                continue

            parent = new_parent

        # Register all tags currently taken by atoms in parent
        temp_taken_tags_uncontested = get_tags_mol(parent)

        # Newly introduced atoms by one of the uncontested reactions need a unique tag
        for atom in parent.GetAtoms():
            if atom.GetIsotope() == 0:  # newly introduced atom has tag 0
                # Loop until we find a tag that is not already taken
                tag = 1
                while tag in (temp_taken_tags_uncontested | original_taken_tags | temp_taken_tags):
                    tag += 1
                atom.SetIsotope(tag)
                temp_taken_tags_uncontested.add(tag)

        unmasked_parent = Mol(parent)  # keep original parent for later
        results = rl.apply(parent, msk)  # apply reaction rule

        try:
            if len(results) == 0:
                raise ValueError(f"No products from uncontested reaction {rl.name}")

            if len(results) > 1:
                raise ValueError(f"More than one product from uncontested reaction {rl.name}")

            result = results[0]
            applied_reactions.append((rl, match_set))  # keep track of successfully applied reactions

            # Reset atom tags in products for atoms not in original reactant
            for product in result:
                for atom in product.GetAtoms():
                    if atom.GetIsotope() not in original_taken_tags and atom.GetIsotope() != 0:
                        atom.SetIsotope(0)
                products.append(product)

        except Exception:
            # Start function again with the next uncontested reaction
            for atom in parent.GetAtoms():
                if atom.GetIsotope() not in original_taken_tags and atom.GetIsotope() != 0:
                    atom.SetIsotope(0)
            products.append(unmasked_parent)
            failed_combos.add((rl.id, frozenset(match_set)))

    for product in products:
        # Any tag in product that is not in parent should be 0; otherwise we run into issues with
        # the set cover algorithm
        for atom in product.GetAtoms():
            if atom.GetIsotope() not in tags_in_parent and atom.GetIsotope() != 0:
                atom.SetIsotope(0)

    return products, applied_reactions, failed_combos


@dataclass(frozen=True)
class MatchingRule:
    """
    Represents a molecular matching rule defined by a SMILES pattern.

    :var name: str: name of the matching rule
    :var smiles: str: SMILES pattern defining the motif
    :var props: dict[str, Any]: additional properties associated with the rule
    :var terminal: bool: whether this rule is terminal (i.e., should not be expanded further)
    :var family_tokens: tuple[str, ...]: tokens representing the family of the matching rule
    :var ancestor_tokens: tuple[tuple[str, ...]]: tokens representing the ancestors of the matching rule
    """

    name: str
    smiles: str
    props: dict[str, Any]
    terminal: bool = True

    family_tokens: set[str] = field(default_factory=set)
    ancestor_tokens: list[str] = field(default_factory=list)

    mol: Mol = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """
        Initialize the Mol from the SMILES pattern.
        """
        mol = smiles_to_mol(self.smiles)
        object.__setattr__(self, "mol", mol)

    @property
    def id(self) -> str:
        """
        Unique identifier for the matching rule based on its SMILES pattern.

        :return: str: unique identifier
        """
        return hashlib.sha256(self.smiles.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the MatchingRule to a dictionary.

        :return: dictionary representation of the MatchingRule
        """
        return {
            "name": self.name,
            "smiles": self.smiles,
            "props": self.props,
            "terminal": self.terminal,
            "family_tokens": list(self.family_tokens),
            "ancestor_tokens": list(self.ancestor_tokens),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "MatchingRule":
        """
        Create a MatchingRule instance from a dictionary.

        :param data: dict[str, Any]: dictionary containing rule data
        :return: MatchingRule: the created MatchingRule instance
        """
        matching_rule = cls(
            name=data["name"],
            smiles=data["smiles"],
            props=data.get("props", {}),
            terminal=data.get("terminal", True),
            family_tokens=set(data.get("family_tokens", [])),
            ancestor_tokens=list(data.get("ancestor_tokens", [])),
        )
        return matching_rule

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

    :var match_stereochemistry: bool: whether to consider stereochemistry in matching rules
    :var reaction_rules: list[ReactionRule]: list of reaction rules
    :var matching_rules: list[MatchingRule]: list of matching rules
    """

    match_stereochemistry: bool
    reaction_rules: list[ReactionRule]
    matching_rules: list[MatchingRule]

    def __str__(self) -> str:
        """
        String representation of the RuleSet.

        :return: str: string representation
        """
        return f"RuleSet({len(self.reaction_rules)} reaction rules, {len(self.matching_rules)} matching rules, match_stereochemistry={self.match_stereochemistry})"

    @classmethod
    def load_default(cls, match_stereochemistry: bool = False) -> "RuleSet":
        """
        Load the default set of reaction and matching rules.

        :return: RuleSet: the default rule set
        """
        path_reaction_rules = Path(files(retromol.data).joinpath("rxn.yml"))
        path_matching_rules_other = Path(files(retromol.data).joinpath("mxn_other.yml"))

        if match_stereochemistry:
            path_matching_rules_polyketide = Path(files(retromol.data).joinpath("mxn_pks_chiral.yml"))
        else:
            path_matching_rules_polyketide = Path(files(retromol.data).joinpath("mxn_pks.yml"))

        with open(path_reaction_rules, "r") as fo:
            reaction_rules_data = yaml.safe_load(fo)

        with open(path_matching_rules_other, "r") as fo:
            matching_rules_other_data = yaml.safe_load(fo)

        with open(path_matching_rules_polyketide, "r") as fo:
            matching_rules_polyketide_data = yaml.safe_load(fo)

        matching_rules_data = matching_rules_other_data + matching_rules_polyketide_data

        reaction_rules = [ReactionRule.from_dict(d) for d in reaction_rules_data]
        matching_rules = [MatchingRule.from_dict(d) for d in matching_rules_data]

        return RuleSet(match_stereochemistry, reaction_rules, matching_rules)
