"""Module defining reaction and matching rules."""

import logging
import itertools
import hashlib
import json
from collections import Counter
from dataclasses import dataclass, field, replace
from functools import cached_property
from importlib.resources import files
from pathlib import Path
from typing import Any, Literal

import yaml
from rdkit.Chem.rdchem import Mol
from rdkit.Chem.rdChemReactions import ChemicalReaction

import retromol.data
from retromol.chem.mol import (
    mol_to_smiles,
    smiles_to_mol,
    smarts_to_mol,
    count_fragments,
    sanitize_mol,
    reassign_stereochemistry,
)
from retromol.chem.reaction import smarts_to_reaction, reactive_template_atoms
from retromol.chem.stereo import BondStereoRecord, capture_double_bond_stereo, restore_double_bond_stereo
from retromol.chem.tagging import get_tags_mol
from retromol.chem.masking import is_masked_preserved


log = logging.getLogger(__name__)

SequenceRole = Literal["backbone", "tailoring"]


@dataclass(frozen=True)
class ReactionRule:
    """
    Represents a chemical reaction rule defined by a SMARTS pattern.

    :cvar name: str: name of the reaction rule
    :cvar smarts: str: SMARTS pattern defining the reaction
    :cvar props: Descriptive metadata, such as notes and references.
    :cvar allowed_in_bulk: bool: whether this rule is allowed to be applied in bulk preprocessing
    :cvar product_max_matches: SMARTS patterns mapped to maximum occurrences in each product.
    """

    name: str
    smarts: str
    props: dict[str, Any]
    allowed_in_bulk: bool = False
    product_max_matches: dict[str, int] = field(default_factory=dict)

    rxn: ChemicalReaction = field(init=False, repr=False)

    def __post_init__(self) -> None:
        """
        Initialize the ChemicalReaction from the SMARTS pattern.
        """
        rxn = smarts_to_reaction(self.smarts)
        object.__setattr__(self, "rxn", rxn)
        # Compile declarative guards when loading, so malformed rules fail early.
        _ = self.product_match_limits

    @cached_property
    def product_match_limits(self) -> tuple[tuple[Mol, int], ...]:
        """Limits on SMARTS occurrences in each product, supplied by the rule."""
        limits = self.product_max_matches
        if not isinstance(limits, dict) or any(
            not isinstance(pattern, str) or not pattern
            or type(limit) is not int or limit < 0
            for pattern, limit in limits.items()
        ):
            raise ValueError(f"{self.name}: product_max_matches must map SMARTS strings to nonnegative integers")
        return tuple((smarts_to_mol(pattern), limit) for pattern, limit in limits.items())

    @property
    def id(self) -> str:
        """
        Unique identifier for the reaction rule based on its SMARTS pattern.

        :return: Unique identifier.
        """
        return hashlib.sha256(self.smarts.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "smarts": self.smarts,
            "props": self.props,
            "allowed_in_bulk": self.allowed_in_bulk,
            "product_max_matches": self.product_max_matches,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "ReactionRule":
        """
        Create a ReactionRule instance from a dictionary.

        :param data: Dictionary containing rule data.
        :return: The created ReactionRule instance.
        """
        # Read older rule files and saved results; serialize only the top-level field.
        props = dict(data.get("props", {}))
        legacy_limits = props.pop("product_max_matches", {})
        reaction_rule = cls(
            name=data["name"],
            smarts=data["smarts"],
            props=props,
            allowed_in_bulk=data.get("allowed_in_bulk", False),
            product_max_matches=data.get("product_max_matches", legacy_limits),
        )
        return reaction_rule
    
    def apply(
        self,
        reactant: Mol,
        mask_tags: set[int] | None = None,
        stereo_registry: dict[frozenset[int], BondStereoRecord] | None = None,
    ) -> list[list[Mol]]:
        """
        Apply the reaction to the given reactant molecule, optionally enforcing a mask on atom tags.

        :param reactant: The reactant molecule.
        :param mask_tags: Set of atom tags (isotope-based tags) that are allowed to change.
        :param stereo_registry: Registry of double bond stereo recorded from the original
            input molecule (see `capture_double_bond_stereo`), used to restore E/Z stereo
            on products that lost it during this reaction. Skipped when None/empty.
        :return: List of unique product tuples (each tuple as a list[Mol]).
        """
        log.debug(f"Applying reaction rule '{self.name}'")

        results = self.rxn.RunReactants([reactant])
        if not results:
            log.debug("No products generated for reactant")
            return []
        log.debug(f"Generated {len(results)} raw product tuple(s)")
        
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
                    log.debug("Product has multiple components, skipping")
                    ok_tuple = False
                    break
                
                # Sanitize in place
                if not sanitize_mol(prod, fix_hydrogens=True):
                    log.debug("Product sanitization failed, skipping")
                    ok_tuple = False
                    break

                if any(len(prod.GetSubstructMatches(pattern, maxMatches=limit + 1)) > limit
                       for pattern, limit in self.product_match_limits):
                    ok_tuple = False
                    break

                # Reassign stereo on the sanitized product
                prod = reassign_stereochemistry(prod)

                # Restore E/Z stereo the reaction couldn't reconstruct itself, where the
                # original input molecule had it defined for this same bond
                if stereo_registry is not None:
                    if stereo_registry:
                        prod = restore_double_bond_stereo(prod, stereo_registry)

                    # A bond can also become stereo-defined by this reaction (e.g. a
                    # ring-opening rule whose product SMARTS hardcodes `/`/`\` geometry
                    # on a bond that was inside a ring, and therefore not stereo-defined,
                    # in the reactant). Record any such new bonds so later reactions on
                    # this product's descendants, which may not hardcode stereo
                    # themselves, can still restore it via the registry.
                    for key, record in capture_double_bond_stereo(prod).items():
                        stereo_registry.setdefault(key, record)

                products.append(prod)
                atom_tag_sets.append(get_tags_mol(prod))

            if not ok_tuple:
                log.debug("Product tuple failed validation, skipping")
                continue

            # Disallow overlapping tag sets across products
            total_tags = sum(len(s) for s in atom_tag_sets)
            union_tags = len(set().union(*atom_tag_sets)) if atom_tag_sets else 0
            if atom_tag_sets and total_tags != union_tags:
                log.debug("Products share atom tags, skipping")
                continue

            # Mask check
            if mask_tags is not None and not is_masked_preserved(reactant, products, mask_tags):
                log.debug("Products modify tags outside mask, skipping")
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

    :param mol: RDKit molecule.
    :param rules: List of preprocessing rules.
    :param failed_combos: Set of failed combinations to avoid infinite loops.
    :return: Uncontested reactions.
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
    stereo_registry: dict[frozenset[int], BondStereoRecord] | None = None,
) -> tuple[list[Mol], list[tuple[ReactionRule, set[int]]], set[tuple[int, frozenset[int]]]]:
    """
    Apply uncontested reactions in bulk.

    :param parent: RDKit molecule.
    :param uncontested: List of uncontested reactions.
    :param original_taken_tags: List of atom tags from original reactant.
    :param stereo_registry: Registry of double bond stereo recorded from the original
        input molecule, forwarded to `ReactionRule.apply` for stereo restoration.
    :return: List of true products, a list of applied ReactionRules with their masks, and a set of failed combinations.
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
        raise ValueError("Not all atoms have a unique tag before applying uncontested reactions!")

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
        results = rl.apply(parent, msk, stereo_registry)  # apply reaction rule

        try:
            if len(results) == 0:
                raise ValueError(f"No products from uncontested reaction {rl.name}!")

            if len(results) > 1:
                raise ValueError(f"More than one product from uncontested reaction {rl.name}!")

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
class MatchingCase:
    """Conditional name/role overrides evaluated after a structural match.

    Omitted or null overrides inherit the matching rule's defaults.
    """

    requires_smarts: str
    name: str | None = None
    sequence_role: SequenceRole | None = None
    _pattern: Mol = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        if not isinstance(self.requires_smarts, str) or not self.requires_smarts.strip():
            raise ValueError("Matching case requires_smarts must be a nonempty SMARTS string")
        try:
            pattern = smarts_to_mol(self.requires_smarts)
        except ValueError as exc:
            raise ValueError(f"Invalid matching case requires_smarts: {self.requires_smarts!r}") from exc
        object.__setattr__(self, "_pattern", pattern)
        if self.name is not None and (not isinstance(self.name, str) or not self.name.strip()):
            raise ValueError("Matching case name must be a nonempty string or null")
        if self.sequence_role not in (None, "backbone", "tailoring"):
            raise ValueError("Matching case sequence_role must be 'backbone', 'tailoring', or null")

    def matches(self, mol: Mol, match_stereochemistry: bool = False) -> bool:
        return mol.HasSubstructMatch(self._pattern, useChirality=match_stereochemistry)

    def to_dict(self) -> dict[str, Any]:
        data = {"requires_smarts": self.requires_smarts}
        if self.name is not None:
            data["name"] = self.name
        if self.sequence_role is not None:
            data["sequence_role"] = self.sequence_role
        return data

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "MatchingCase":
        if not isinstance(data, dict) or "requires_smarts" not in data:
            raise ValueError("Each matching case must be a mapping with requires_smarts")
        unknown = data.keys() - {"requires_smarts", "name", "sequence_role"}
        if unknown:
            raise ValueError(f"Unknown matching case fields: {sorted(unknown)}")
        return cls(**data)


@dataclass(frozen=True)
class MatchingRule:
    """
    Represents a molecular matching rule defined by a SMILES pattern.

    :cvar name: Name of the matching rule.
    :cvar smiles: SMILES pattern defining the motif.
    :cvar props: Descriptive metadata, such as notes and references.
    :cvar pseudonyms: Pseudonyms associated with the rule.
    :cvar terminal: Whether this rule is terminal (i.e., should not be expanded further).
    :cvar stereochemistry: Whether to consider stereochemistry in matching.
    :cvar display_smiles: Optional friendlier SMILES to show when depicting this motif (e.g.
        in the GUI's motif hover preview) instead of `smiles`, which is written for
        substructure matching (e.g. carrying a reactive leaving-group placeholder) and
        may not be the clearest depiction of the motif on its own.
    :cvar preserve_across_routes: Whether a terminal match also protects its atom
        region from subdivision in alternative derivations. Defaults to True.
    :cvar sequence_role: Optional role in a linear assembly readout. None uses
        the bond-element heuristic; this setting does not affect parsing.
    :cvar requires_smarts: Optional additional substructure query evaluated on
        the tagged fragment after the whole-fragment SMILES match. Isotope zero
        can distinguish reaction-added atoms from original input atoms.
    :cvar cases: Conditional name/role overrides, evaluated only after this rule
        matches. The first matching case wins; otherwise the defaults apply.
    """

    name: str
    smiles: str
    props: dict[str, Any]
    pseudonyms: list[str]
    terminal: bool = True
    stereochemistry: bool = False
    display_smiles: str | None = None
    preserve_across_routes: bool = True
    sequence_role: SequenceRole | None = None
    requires_smarts: str | None = None
    cases: list[MatchingCase] = field(default_factory=list)

    mol: Mol = field(init=False, repr=False)
    _required_pattern: Mol | None = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        """
        Initialize the Mol from the SMILES pattern.
        """
        mol = smiles_to_mol(self.smiles)
        object.__setattr__(self, "mol", mol)
        if type(self.preserve_across_routes) is not bool:
            raise ValueError(f"{self.name}: preserve_across_routes must be a boolean")
        if self.sequence_role not in (None, "backbone", "tailoring"):
            raise ValueError(f"{self.name}: sequence_role must be 'backbone', 'tailoring', or null")
        pattern = None
        if self.requires_smarts is not None:
            if not isinstance(self.requires_smarts, str) or not self.requires_smarts.strip():
                raise ValueError(f"{self.name}: requires_smarts must be a nonempty SMARTS string or null")
            try:
                pattern = smarts_to_mol(self.requires_smarts)
            except ValueError as exc:
                raise ValueError(f"{self.name}: invalid requires_smarts: {self.requires_smarts!r}") from exc
        object.__setattr__(self, "_required_pattern", pattern)
        if not isinstance(self.cases, list) or any(not isinstance(case, MatchingCase) for case in self.cases):
            raise ValueError(f"{self.name}: cases must be a list of MatchingCase objects")

    @property
    def id(self) -> str:
        """
        Identifier for the SMILES and optional SMARTS condition.

        Rules without a condition retain their historical SMILES-only IDs.

        :return: Unique identifier.
        """
        definition = self.smiles if self.requires_smarts is None else json.dumps(
            {"smiles": self.smiles, "requires_smarts": self.requires_smarts},
            sort_keys=True, separators=(",", ":"),
        )
        return hashlib.sha256(definition.encode("utf-8")).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        """
        Serialize the MatchingRule to a dictionary.

        :return: Dictionary representation of the MatchingRule.
        """
        return {
            "name": self.name,
            "smiles": self.smiles,
            "props": self.props,
            "pseudonyms": self.pseudonyms,
            "terminal": self.terminal,
            "stereochemistry": self.stereochemistry,
            "display_smiles": self.display_smiles,
            "preserve_across_routes": self.preserve_across_routes,
            "sequence_role": self.sequence_role,
            "requires_smarts": self.requires_smarts,
            "cases": [case.to_dict() for case in self.cases],
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "MatchingRule":
        """
        Create a MatchingRule instance from a dictionary.

        :param data: Dictionary containing rule data.
        :return: The created MatchingRule instance.
        """
        # Keep old YAML/JSON readable without retaining behavior settings in metadata.
        props = dict(data.get("props", {}))
        legacy_preserve = props.pop("preserve_across_routes", True)
        cases = data.get("cases", [])
        if not isinstance(cases, list):
            raise ValueError(f"{data['name']}: cases must be a list")
        matching_rule = cls(
            name=data["name"],
            smiles=data["smiles"],
            props=props,
            pseudonyms=data.get("pseudonyms", []),
            terminal=data.get("terminal", True),
            stereochemistry=data.get("stereochemistry", False),
            display_smiles=data.get("display_smiles"),
            preserve_across_routes=data.get("preserve_across_routes", legacy_preserve),
            sequence_role=data.get("sequence_role"),
            requires_smarts=data.get("requires_smarts"),
            cases=[MatchingCase.from_dict(case) for case in cases],
        )
        return matching_rule

    def is_match(self, mol: Mol, match_stereochemistry: bool = False) -> bool:
        """
        Check if the given molecule matches this rule.

        :param mol: Molecule to check.
        :param match_stereochemistry: Whether to consider stereochemistry in matching.
        :return: True if the molecule matches the rule, False otherwise.
        """
        has_substruct_match = mol.HasSubstructMatch(self.mol, useChirality=match_stereochemistry)
        has_equal_num_atoms = mol.GetNumAtoms() == self.mol.GetNumAtoms()
        has_equal_num_bonds = mol.GetNumBonds() == self.mol.GetNumBonds()
        
        if has_substruct_match and has_equal_num_atoms and has_equal_num_bonds:
            return self._required_pattern is None or mol.HasSubstructMatch(
                self._required_pattern, useChirality=match_stereochemistry,
            )

        return False

    def _resolved(self, case: MatchingCase | None = None) -> "MatchingRule":
        """Snapshot the effective identity without modifying the source rule."""
        if not self.cases:
            return self
        return replace(
            self,
            name=case.name if case is not None and case.name is not None else self.name,
            sequence_role=(case.sequence_role if case is not None and case.sequence_role is not None
                           else self.sequence_role),
            cases=[],
        )

    def match(self, mol: Mol, match_stereochemistry: bool = False) -> "MatchingRule | None":
        """Match the structure, then resolve its name and role using cases."""
        if not self.is_match(mol, match_stereochemistry):
            return None
        for case in self.cases:
            if case.matches(mol, match_stereochemistry):
                return self._resolved(case)
        return self._resolved()

    def variants(self) -> list["MatchingRule"]:
        """Possible identities for name lookup/depiction, not parser matching."""
        return [self._resolved(), *(self._resolved(case) for case in self.cases)]


@dataclass(frozen=True)
class RuleSet:
    """
    Represents a set of reaction and matching rules.

    :var match_stereochemistry: Whether to consider stereochemistry in matching rules.
    :var reaction_rules: List of reaction rules.
    :var matching_rules: List of matching rules.
    """

    match_stereochemistry: bool
    reaction_rules: list[ReactionRule]
    matching_rules: list[MatchingRule]

    @cached_property
    def reaction_rules_hash(self) -> str:
        """SHA-256 of ordered, canonical rule definitions (not YAML formatting).

        Treat a RuleSet and its nested lists/props as immutable once used.
        """
        return self._hash_rules(self.reaction_rules)

    @cached_property
    def matching_rules_hash(self) -> str:
        """SHA-256 of ordered matching definitions, including names and flags."""
        return self._hash_rules(self.matching_rules)

    @cached_property
    def matching_rule_variants(self) -> list[MatchingRule]:
        """Default and case identities for vocabularies and name-based lookup.

        Parsing always uses matching_rules, so each structure is matched once
        before evaluating its cases. Variants may share a structural rule ID.
        """
        return [variant for rule in self.matching_rules for variant in rule.variants()]

    @staticmethod
    def _hash_rules(rules: list) -> str:
        content = json.dumps(
            [rule.to_dict() for rule in rules], sort_keys=True,
            separators=(",", ":"), ensure_ascii=False, allow_nan=False,
        )
        return hashlib.sha256(content.encode("utf-8")).hexdigest()

    def __str__(self) -> str:
        """
        String representation of the RuleSet.

        :return: String representation.
        """
        return f"RuleSet({len(self.reaction_rules)} reaction rules, {len(self.matching_rules)} matching rules, match_stereochemistry={self.match_stereochemistry})"

    def find_structural_matches(self, mol: Mol, match_stereochemistry: bool = False) -> list["MatchingRule"]:
        """
        Find every matching rule and resolve its cases for this tagged molecule.

        :param mol: Molecule to match against every rule in this ruleset.
        :param match_stereochemistry: Whether to consider stereochemistry in matching.
        :return: Every MatchingRule whose motif has the same structure as `mol`, in ruleset order (possibly empty).
        """
        matches = []
        for rule in self.matching_rules:
            resolved = rule.match(mol, match_stereochemistry)
            if resolved is not None:
                matches.append(resolved)
        return matches

    @classmethod
    def load_from_files(
        cls,
        reaction_rules_path: str | Path | None = None,
        matching_rules_path: str | Path | None = None,
        match_stereochemistry: bool = False,
    ) -> "RuleSet":
        """
        Load a set of reaction and matching rules from a file.

        :param reaction_rules_path: Path to file containing reaction_rules.
        :param matching_rules_path: Path to file containing matching rules.
        :param match_stereochemistry: Whether to consider stereochemistry in matching rules.
        :return: A set of reaction rules and matching rules.
        """
        path_reaction_rules = Path(reaction_rules_path) if reaction_rules_path else Path(files(retromol.data).joinpath("rxn.yml"))
        path_matching_rules = Path(matching_rules_path) if matching_rules_path else Path(files(retromol.data).joinpath("mxn.yml"))

        with open(path_reaction_rules, "r") as fo:
            reaction_rules_data = yaml.safe_load(fo)

        with open(path_matching_rules, "r") as fo:
            matching_rules_data = yaml.safe_load(fo)

        reaction_rules = [ReactionRule.from_dict(d) for d in reaction_rules_data]
        matching_rules = [MatchingRule.from_dict(d) for d in matching_rules_data]

        matching_rules.sort(key=lambda rule: not rule.stereochemistry)

        return RuleSet(match_stereochemistry, reaction_rules, matching_rules)

    @classmethod
    def load_default(cls, match_stereochemistry: bool = False) -> "RuleSet":
        """
        Load the default RuleSet from package resources.

        :param match_stereochemistry: Whether to consider stereochemistry in matching rules.
        :return: The loaded RuleSet instance.
        """
        return cls.load_from_files(match_stereochemistry=match_stereochemistry)
