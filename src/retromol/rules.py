"""This module contains functions for parsing RetroMol rules."""

import logging
from collections import defaultdict
from copy import deepcopy
from dataclasses import dataclass
from importlib.resources import files
from typing import Any

import yaml
from rdkit.Chem.rdchem import PeriodicTable
from rdkit.Chem.rdMolDescriptors import CalcNumRings
from rdkit.Chem.rdmolops import (
    AssignAtomChiralTagsFromStructure,
    AssignStereochemistry,
    SetBondStereoFromDirections,
)
from tqdm import tqdm

import retromol.data
from retromol.chem import (
    ChemicalReaction,
    Mol,
    count_fragments,
    get_default_valence,
    get_periodic_table,
    get_tags_mol,
    mol_to_inchikey,
    mol_to_smiles,
    sanitize_mol,
    smarts_to_mol,
    smarts_to_reaction,
    smiles_to_mol,
    stereo_summary,
)
from retromol.config import LOGGER_NAME
from retromol.helpers import sha256_hex


def check_tags_are_nonzero(mol: Mol) -> None:
    """
    Check if all atom tags are nonzero.

    :param rea: molecule

    :raises ValueError: if any atom tag is 0
    """
    curr_tags = [a.GetIsotope() for a in mol.GetAtoms()]

    if any([tag == 0 for tag in curr_tags]):
        raise ValueError("molecule contains atom tag 0")


def apply_mask(mol: Mol, msk: set[int]) -> dict[int, int]:
    """
    Set atom numbers of atoms not in mask to 0 based on atom tags.

    :param mol: molecule
    :param msk: mask of atom tags
    :return: mapping of atom tags to atomic numbers
    .. note:: this function modifies the reactant molecule in place
    """
    # Check if all atom tags are nonzero
    check_tags_are_nonzero(mol)

    # Get original tag to atomic num mapping
    tag_to_anr = {a.GetIsotope(): a.GetAtomicNum() for a in mol.GetAtoms()}

    # Now we can apply the mask, set atomic num to 0 for atoms not in mask
    for atom in mol.GetAtoms():
        if atom.GetIsotope() not in msk:
            atom.SetAtomicNum(0)

    return tag_to_anr


def check_atomic_nums_are_nonzero(mol: Mol) -> None:
    """
    Check if all atomic numbers are nonzero.

    :param rea: molecule
    :raises ValueError: if any atomic number is 0
    """
    curr_anrs = [a.GetAtomicNum() for a in mol.GetAtoms()]

    if any([anr == 0 for anr in curr_anrs]):
        raise ValueError("molecule contains atomic number 0")


def reset_atomic_nums(mol: Mol, tag_to_anr: dict[int, int]) -> None:
    """
    Reset atomic numbers of atoms based on atom tags.

    :param mol: molecule
    :param tag_to_anr: mapping of atom tags to atomic numbers
    .. note:: this function modifies the reactant molecule in place
    """
    for atom in mol.GetAtoms():
        if atom.GetIsotope() != 0:  # newly added atoms in reaction have isotope 0
            original_anr = tag_to_anr.get(atom.GetIsotope(), None)
            if original_anr is None:
                raise ValueError(f"no atomic num found for atom tag {atom.GetIsotope()}")
            atom.SetAtomicNum(original_anr)

    # Make sure that no atomic num is 0
    check_atomic_nums_are_nonzero(mol)


def correct_hydrogens(mol: Mol) -> None:
    """
    Correct explicit hydrogens on atoms based on valence rules.

    :param mol: molecule
    .. note:: this function modifies the reactant molecule in place
    """
    for atom in mol.GetAtoms():
        # Skip aromatic and charged atoms
        if atom.GetIsAromatic() or atom.GetFormalCharge() != 0:
            continue

        # Skip phosphorus and sulfur (can have expanded valence)
        if atom.GetAtomicNum() in (15, 16):
            continue

        # Check if atom complies with valence rules, otherwise adjust explicit Hs
        valence_bonds = int(sum([bond.GetValenceContrib(atom) for bond in atom.GetBonds()]))
        default_valence = get_default_valence(atom.GetAtomicNum())
        num_hs = atom.GetNumExplicitHs()

        if default_valence - valence_bonds < num_hs:
            new_valence = default_valence - valence_bonds

            if new_valence < 0:
                raise ValueError("new atom valence is negative")

            atom.SetNumExplicitHs(new_valence)


def _collect_map_to_atomicnums(mols: list[Mol]) -> dict[int, list[int]]:
    """
    Return map_number -> list of atomic numbers (with multiplicity).

    :param mols: list of molecules
    :return: mapping of map numbers to list of atomic numbers
    """
    out: dict[int, list[int]] = defaultdict(list)
    for mol in mols:
        for atom in mol.GetAtoms():
            m = atom.GetAtomMapNum()
            if m > 0:
                out[m].append(atom.GetAtomicNum())
    return out


def summarize_atomicnums(nums: list[int]) -> str:
    """
    Compact helper for error messages, e.g., [6,6,8] -> Cx2,O.

    :param nums: list of atomic numbers
    :return: compact string representation
    """
    if not nums:
        return "[]"
    syms = [PeriodicTable.GetElementSymbol(get_periodic_table(), z) for z in nums]
    # Build counts
    from collections import Counter

    c = Counter(syms)
    parts = [f"{el}x{cnt}" if cnt > 1 else el for el, cnt in sorted(c.items())]
    return "[" + ",".join(parts) + "]"


@dataclass(frozen=True)
class _CompiledConds:
    """
    Compiled conditions for reaction rules.

    :param requires_any: list of required substructures (any)
    :param requires_all: list of required substructures (all)
    :param forbids_any: list of forbidden substructures (any)
    :param min_counts: list of (substructure, min count) tuples
    :param max_counts: list of (substructure, max count) tuples
    :param ring_min: minimum number of rings
    :param ring_max: maximum number of rings
    :param atom_min: minimum number of atoms
    :param atom_max: maximum number of atoms
    :param charge_min: minimum total charge
    :param charge_max: maximum total charge
    :param has_metal: whether the molecule must contain a metal
    :param is_macrocycle: whether the molecule must be a macrocycle
    """

    requires_any: list[Mol]
    requires_all: list[Mol]
    forbids_any: list[Mol]
    min_counts: list[tuple[Mol, int]]
    max_counts: list[tuple[Mol, int]]
    ring_min: int | None
    ring_max: int | None
    atom_min: int | None
    atom_max: int | None
    charge_min: int | None
    charge_max: int | None
    has_metal: bool | None
    is_macrocycle: bool | None


def _compile_smarts_list(lst: list[str] | None) -> list[Mol]:
    """
    Compile a list of SMARTS strings into RDKit Mol objects.

    :param lst: list of SMARTS strings
    :return: list of RDKit Mol objects
    """
    if not lst:
        return []
    out: list[Mol] = []
    for s in lst:
        m = smarts_to_mol(s)
        out.append(m)
    return out


def _compile_counts(d: dict[str, int] | None) -> list[tuple[Mol, int]]:
    """
    Compile a dictionary of SMARTS strings to counts into a list of (Mol, count) tuples.

    :param d: dictionary of SMARTS strings to counts
    :return: list of (Mol, count) tuples
    """
    if not d:
        return []
    return [(smarts_to_mol(k), v) for k, v in d.items()]


# Set of atomic numbers considered as metals
_METALS = {
    3,
    11,
    19,
    37,
    55,
    87,
    4,
    12,
    20,
    38,
    56,
    88,
    *range(21, 31),
    *range(39, 49),
    *range(72, 81),
    *range(57, 72),
    *range(89, 104),
    13,
    31,
    49,
    50,
    81,
    82,
    83,
}


def _has_metal(mol: Mol) -> bool:
    """
    Check if the molecule contains any metal atoms.

    :param mol: molecule
    :return: True if the molecule contains metal atoms, otherwise False
    """
    return any(a.GetAtomicNum() in _METALS for a in mol.GetAtoms())


def _has_macrocycle(mol: Mol, min_size: int = 12) -> bool:
    """
    Check if the molecule contains a macrocycle (ring of at least min_size).

    :param mol: molecule
    :param min_size: minimum size of the ring to be considered a macrocycle
    :return: True if the molecule contains a macrocycle, otherwise False
    """
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) >= min_size:
            return True
    return False


def _passes_global(mol: Mol, C: _CompiledConds) -> bool:
    """
    Check if a molecule passes the global conditions.

    :param mol: molecule
    :param C: compiled conditions
    :return: True if the molecule passes the conditions, otherwise False
    """
    # SMARTS presence/absence
    if C.requires_any and not any(mol.HasSubstructMatch(q) for q in C.requires_any):
        return False
    if any(not mol.HasSubstructMatch(q) for q in C.requires_all):
        return False
    if any(mol.HasSubstructMatch(q) for q in C.forbids_any):
        return False

    # Count thresholds
    for q, n in C.min_counts:
        if len(mol.GetSubstructMatches(q)) < n:
            return False
    for q, n in C.max_counts:
        if len(mol.GetSubstructMatches(q)) > n:
            return False

    # Simple numeric props
    n_ring = CalcNumRings(mol)
    if C.ring_min is not None and n_ring < C.ring_min:
        return False
    if C.ring_max is not None and n_ring > C.ring_max:
        return False

    n_atoms = mol.GetNumAtoms()
    if C.atom_min is not None and n_atoms < C.atom_min:
        return False
    if C.atom_max is not None and n_atoms > C.atom_max:
        return False

    charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    if C.charge_min is not None and charge < C.charge_min:
        return False
    if C.charge_max is not None and charge > C.charge_max:
        return False

    if C.has_metal is not None and _has_metal(mol) != C.has_metal:
        return False
    if C.is_macrocycle is not None and _has_macrocycle(mol) != C.is_macrocycle:
        return False

    return True


@dataclass(frozen=True)
class ReactionRule:
    """
    Preprocessing rule for a chemical reaction.

    :param id: internal identifier
    :param rid: rule identifier
    :param rxn: RDKit ChemicalReaction (optional if `smarts` is provided)
    :param smarts: reaction SMARTS string
    :param groups: groups the reaction rule belongs to
    :param props: properties
    """

    id: int
    rid: str
    rxn: ChemicalReaction | None
    smarts: str
    groups: list[str]
    props: dict[str, Any]

    def __post_init__(self) -> None:
        # Allow constructing rxn from smarts in a frozen class
        rxn = self.rxn
        if rxn is None:
            if not self.smarts:
                raise ValueError(f"[{self.rid}] Either rxn or smarts must be provided.")
            try:
                rxn = smarts_to_reaction(self.smarts, use_smiles=False)
            except Exception as e:
                raise ValueError(f"[{self.rid}] Invalid reaction SMARTS: {e}") from e
            object.__setattr__(self, "rxn", rxn)

        # Basic sanity: at least one reactant and one product
        if rxn.GetNumReactantTemplates() == 0 or rxn.GetNumProductTemplates() == 0:
            raise ValueError(f"[{self.rid}] Reaction must have >=1 reactant and >=1 product.")

        # Collect mapped atoms (molAtomMapNumber) from both sides
        reactant_maps = _collect_map_to_atomicnums(list(rxn.GetReactants()))
        product_maps = _collect_map_to_atomicnums(list(rxn.GetProducts()))

        # Keys (map numbers) must match exactly
        maps_react = set(reactant_maps.keys())
        maps_prod = set(product_maps.keys())
        if maps_react != maps_prod:
            missing_in_prod = sorted(maps_react - maps_prod)
            missing_in_reac = sorted(maps_prod - maps_react)
            msgs: list[str] = []
            if missing_in_prod:
                msgs.append(f"map nums only on reactant side: {missing_in_prod}")
            if missing_in_reac:
                msgs.append(f"map nums only on product side: {missing_in_reac}")
            raise ValueError(f"[{self.rid}] Mapped atoms don't add up: " + "; ".join(msgs))

        # For every map number, multiplicity and atomic numbers must match
        errors: list[str] = []
        for m in sorted(maps_react):
            r_list = sorted(reactant_maps[m])
            p_list = sorted(product_maps[m])
            if r_list != p_list:
                # This captures both multiplicity and element changes.
                errors.append(
                    f"map {m}: reactants {summarize_atomicnums(r_list)} != products {summarize_atomicnums(p_list)}"
                )

        if errors:
            raise ValueError(
                f"[{self.rid}] Mapped atom consistency failed (no element changes allowed): " + "; ".join(errors)
            )

        # Read and compile conditions
        reactant_conds = self.props.get("conditions", {}).get("reactant", {})
        product_conds = self.props.get("conditions", {}).get("product", {})

        object.__setattr__(
            self,
            "_reactant_conds",
            _CompiledConds(
                requires_any=_compile_smarts_list(reactant_conds.get("requires_any")),
                requires_all=_compile_smarts_list(reactant_conds.get("requires_all")),
                forbids_any=_compile_smarts_list(reactant_conds.get("forbids_any")),
                min_counts=_compile_counts(reactant_conds.get("min_counts")),
                max_counts=_compile_counts(reactant_conds.get("max_counts")),
                ring_min=(reactant_conds.get("ring_count") or {}).get("min"),
                ring_max=(reactant_conds.get("ring_count") or {}).get("max"),
                atom_min=(reactant_conds.get("atom_count") or {}).get("min"),
                atom_max=(reactant_conds.get("atom_count") or {}).get("max"),
                charge_min=(reactant_conds.get("total_charge") or {}).get("min"),
                charge_max=(reactant_conds.get("total_charge") or {}).get("max"),
                has_metal=reactant_conds.get("custom_props", {}).get("has_metal"),
                is_macrocycle=reactant_conds.get("custom_props", {}).get("is_macrocycle"),
            ),
        )
        object.__setattr__(
            self,
            "_product_conds",
            _CompiledConds(
                requires_any=_compile_smarts_list(product_conds.get("requires_any")),
                requires_all=_compile_smarts_list(product_conds.get("requires_all")),
                forbids_any=_compile_smarts_list(product_conds.get("forbids_any")),
                min_counts=_compile_counts(product_conds.get("min_counts")),
                max_counts=_compile_counts(product_conds.get("max_counts")),
                ring_min=(product_conds.get("ring_count") or {}).get("min"),
                ring_max=(product_conds.get("ring_count") or {}).get("max"),
                atom_min=(product_conds.get("atom_count") or {}).get("min"),
                atom_max=(product_conds.get("atom_count") or {}).get("max"),
                charge_min=(product_conds.get("total_charge") or {}).get("min"),
                charge_max=(product_conds.get("total_charge") or {}).get("max"),
                has_metal=product_conds.get("custom_props", {}).get("has_metal"),
                is_macrocycle=product_conds.get("custom_props", {}).get("is_macrocycle"),
            ),
        )

    def to_json_serializable_dict(self) -> dict[str, Any]:
        """
        Convert the reaction rule to a JSON serializable dictionary.

        :return: JSON serializable dictionary
        """
        return {
            "rid": self.rid,
            "smarts": self.smarts,
            "groups": self.groups,
            "props": self.props,
        }

    @classmethod
    def from_json_serializable_dict(cls, internal_identifier: int, d: dict[str, Any]) -> "ReactionRule":
        """
        Convert a JSON serializable dictionary to a ReactionRule.

        :param internal_identifier: internal identifier for the reaction rule
        :param d: JSON serializable dictionary
        :return: ReactionRule
        """
        return ReactionRule(
            id=internal_identifier,
            rid=d["rid"],
            rxn=smarts_to_reaction(d["smarts"]),
            smarts=d["smarts"],
            groups=d["groups"],
            props=d.get("props", {}),
        )

    def has_ring_matching_condition(self) -> bool:
        """
        Check if the reaction has a ring matching condition.

        :return: True if the reaction has a ring matching condition, otherwise False
        """
        return any([self.smarts.find(f";{rc}") != -1 for rc in ["R", "!R"]])

    def expected_num_products(self) -> int:
        """
        Get the expected number of products for the reaction.

        :return: expected number of products
        """
        if self.rxn is None:
            raise ValueError("reaction is not initialized")

        return self.rxn.GetNumProductTemplates()

    def __call__(self, rea: Mol, msk: set[int] | None = None) -> list[list[Mol]]:
        """
        Apply the reaction, sanitize, preserve/reassign stereochemistry,
        enforce mask WITHOUT mutating atomic numbers, and dereplicate results
        in a stereo-aware, multiplicity-preserving, order-insensitive way.

        :param rea: reactant molecule
        :param msk: set of atom tags (isotope-based tags) that are allowed to change
        :return: list of unique product tuples (each tuple as a list[Mol])
        """
        logger = logging.getLogger(LOGGER_NAME)

        if self.rxn is None:
            raise ValueError("reaction is not initialized")

        def _prepare_stereo(m: Mol) -> Mol:
            # Reassign stereo cleanly without changing identity
            mm = Mol(m)
            SetBondStereoFromDirections(mm)
            if mm.GetNumConformers() > 0:
                AssignAtomChiralTagsFromStructure(mm, replaceExistingTags=True)
            AssignStereochemistry(mm, cleanIt=True, force=True, flagPossibleStereoCenters=True)
            return mm

        def _single_component(m: Mol) -> bool:
            return count_fragments(m) == 1

        def _sanitize_in_place(m: Mol) -> bool:
            try:
                correct_hydrogens(m)
                sanitize_mol(m)
                return True
            except ValueError:
                return False

        def _tag_to_idx(m: Mol) -> dict[int, int]:
            # "atom tags" live in Isotope; ignore zeros
            d: dict[int, int] = {}
            for a in m.GetAtoms():
                t = a.GetIsotope()
                if t:
                    d[t] = a.GetIdx()
            return d

        def _neighbor_sig(m: Mol, ai: int) -> list[tuple[int, float]]:
            # Neighbor signature by (neighbor tag or neighbor atomicnum if untagged, bond order)
            out: list[tuple[int, float]] = []
            a = m.GetAtomWithIdx(ai)
            for b in a.GetBonds():
                nb = b.GetOtherAtomIdx(ai)
                na = m.GetAtomWithIdx(nb)
                ntag = na.GetIsotope()
                key = ntag if ntag else -na.GetAtomicNum()
                out.append((key, float(b.GetBondTypeAsDouble())))
            out.sort()
            return out

        def _mapped_tags_changed(r: Mol, p: Mol) -> set[int]:
            """
            Heuristic diff between reactant and product by tags.
            A tag is considered 'changed' if:
            - the atom with that tag changes atomic number, OR
            - its neighbor signature (by tagged IDs / atom types + bond order) changes.

            :param r: reactant molecule
            :param p: product molecule
            :return: set of changed atom tags
            """
            changed: set[int] = set()
            rmap = _tag_to_idx(r)
            pmap = _tag_to_idx(p)
            for t in set(rmap).intersection(pmap):
                ra = r.GetAtomWithIdx(rmap[t])
                pa = p.GetAtomWithIdx(pmap[t])
                if ra.GetAtomicNum() != pa.GetAtomicNum():
                    changed.add(t)
                    continue
                if _neighbor_sig(r, rmap[t]) != _neighbor_sig(p, pmap[t]):
                    changed.add(t)
            return changed

        def _preserves_mask(reactant: Mol, products: list[Mol], allowed: set[int]) -> bool:
            """
            Check that only tags in `allowed` are changed across all products.

            :param reactant: reactant molecule
            :param products: list of product molecules
            :param allowed: set of allowed changed tags
            :return: True if only allowed tags are changed, otherwise False
            """
            if not allowed:
                return True
            changed: set[int] = set()
            for pr in products:
                changed |= _mapped_tags_changed(reactant, pr)
            return changed.issubset(allowed)

        def _product_key(m: Mol) -> str:
            """
            Stereo-aware, mapping-/tag-invariant canonical key for a product.
            - Clears molAtomMapNumber props (if present) and strips isotope tags ONLY on the copy
                so symmetric mappings don't duplicate results.
            - Uses isomeric canonical SMILES to preserve R/S and E/Z.

            :param m: product molecule
            :return: product key string
            """
            mm = Mol(m)
            for a in mm.GetAtoms():
                if a.HasProp("molAtomMapNumber"):
                    a.ClearProp("molAtomMapNumber")
                # Strip isotope-based tags for the KEY only:
                if a.GetIsotope():
                    a.SetIsotope(0)
            # Important: isomericSmiles preserves stereo
            return mol_to_smiles(mm, isomeric=True, canonical=True)

        def _result_key(products: list[Mol]) -> tuple[tuple[str, int], ...]:
            """
            Stereo-aware, order-insensitive,

            multiplicity-preserving key for a product tuple.
            :param products: list of product molecules
            :return: result key tuple
            """
            # multiplicity-preserving multiset key (Counter of product keys)
            from collections import Counter

            c = Counter(_product_key(p) for p in products)
            return tuple(sorted(c.items(), key=lambda kv: kv[0]))

        logger.debug(f"({self.rid}) applying reaction rule to... {mol_to_smiles(deepcopy(rea), remove_tags=True)}")

        # Pre-filter on reactant
        if not _passes_global(rea, self._reactant_conds):
            logger.debug(f"({self.rid}) reactant fails global conditions")
            return []

        # Run reaction
        results = self.rxn.RunReactants([rea])
        if not results:
            logger.debug(f"({self.rid}) no valid products found after applying RDKit reaction")
            return []

        # Sanitize, filter invalids, product-conditions
        kept: list[list[Mol]] = []
        for tup in results:
            products: list[Mol] = []

            # Quick shape check + sanitize
            atom_tag_sets: list[set[int]] = []
            ok_tuple = True
            for prod in tup:
                if not _single_component(prod):
                    logger.debug(
                        f"({self.rid}) product not single component: {mol_to_smiles(deepcopy(prod), remove_tags=True)}"
                    )
                    ok_tuple = False
                    break

                if not _sanitize_in_place(prod):
                    logger.debug(
                        f"({self.rid}) product failed sanitization: {mol_to_smiles(deepcopy(prod), remove_tags=True)}"
                    )
                    ok_tuple = False
                    break

                # Reassign stereo on the sanitized product
                prod = _prepare_stereo(prod)

                products.append(prod)
                atom_tag_sets.append(set(get_tags_mol(prod)))

            if not ok_tuple:
                logger.debug(f"({self.rid}) product failed validation")
                continue

            # Disallow overlapping tag sets across products
            total_tags = sum(len(s) for s in atom_tag_sets)
            union_tags = len(set().union(*atom_tag_sets)) if atom_tag_sets else 0
            if atom_tag_sets and total_tags != union_tags:
                logger.debug(f"({self.rid}) products share atom tags: {[mol_to_smiles(p) for p in products]}")
                continue

            # Product-side global conditions
            if not all(_passes_global(p, self._product_conds) for p in products):
                logger.debug(f"({self.rid}) products fail global conditions")
                continue

            # Mask check
            if msk is not None and not _preserves_mask(rea, products, msk):
                logger.debug(f"({self.rid}) products modify tags outside mask")
                # Skip results that modify tags outside the mask
                continue

            kept.append(products)

        if len(kept) <= 1:
            return kept

        # Stereo-aware derep (order-insensitive, multiplicity-preserving)
        seen: dict[tuple[tuple[str, int], ...], int] = {}
        unique: list[list[Mol]] = []
        for res in kept:
            key = _result_key(res)
            if key in seen:
                continue
            seen[key] = 1
            unique.append(res)

        return unique


def DummyReactionRule(rid: str) -> ReactionRule:
    """
    Create a dummy reaction rule for testing purposes.

    :param rid: rule identifier
    :return: ReactionRule
    """
    return ReactionRule(
        id=0,
        rid=rid,
        rxn=smarts_to_reaction("[C:1]>>[C:1]"),
        smarts="[C:1]>>[C:1]",
        groups=[],
        props={},
    )


@dataclass(frozen=True)
class MatchingRule:
    """
    Matching rule for a chemical compound.

    :param id: internal identifier
    :param rid: rule identifier
    :param smiles: SMILES string of the molecule
    :param mol: molecule to match against
    :param groups: groups the matching rule belongs to
    :param props: properties of the rule
    """

    id: int
    rid: str
    smiles: str
    mol: Mol
    groups: list[str]
    props: dict[str, Any]

    def to_json_serializable_dict(self) -> dict[str, Any]:
        """Convert the matching rule to a JSON serializable dictionary.

        :return: JSON serializable dictionary
        """
        return {
            "rid": self.rid,
            "mol": mol_to_smiles(self.mol),
            "groups": self.groups,
            "props": self.props,
        }

    @classmethod
    def from_json_serializable_dict(cls, internal_identifier: int, d: dict[str, Any]) -> "MatchingRule":
        """
        Convert a JSON serializable dictionary to a MatchingRule.

        :param internal_identifier: internal identifier for the matching rule
        :param d: JSON serializable dictionary
        :return: MatchingRule
        """
        return MatchingRule(
            id=internal_identifier,
            rid=d["rid"],
            smiles=d["mol"],
            mol=smiles_to_mol(d["mol"]),
            groups=d.get("groups", []),
            props=d.get("props", {}),
        )

    def is_match(self, mol: Mol, sch: bool = False) -> str | None:
        """
        Check if the molecule matches the rule.

        :param mol: molecule
        :param sch: flag to enable/disable stereochemistry matching
        :return: name of the rule if the molecule matches, otherwise None
        """
        # No stereochemistry matching, check for substructure match
        has_substruct_match = mol.HasSubstructMatch(self.mol, useChirality=sch)
        has_equal_num_atoms = mol.GetNumAtoms() == self.mol.GetNumAtoms()
        has_equal_num_bonds = mol.GetNumBonds() == self.mol.GetNumBonds()
        if has_substruct_match and has_equal_num_atoms and has_equal_num_bonds:
            if sch:
                tag = stereo_summary(mol)
                return f"{self.rid}[{tag}]" if tag != "none" else self.rid
            return self.rid
        else:
            return None


class Rules:
    def __init__(
        self,
        reaction_rules: list[ReactionRule],
        matching_rules: list[MatchingRule],
        sha256_reaction_rules: str | None = None,
        sha256_matching_rules: str | None = None,
    ) -> None:
        """
        Initialize the rules for processing compounds.

        :param reaction_rules: list of reaction rules
        :param matching_rules: list of matching rules
        :param sha256_reaction_rules: SHA256 hash of the reaction rules (optional)
        :param sha256_matching_rules: SHA256 hash of the matching rules (optional)
        """
        self._logger = logging.getLogger(LOGGER_NAME)
        self._reaction_rules = reaction_rules
        self._matching_rules = matching_rules
        self.sha256_reaction_rules = sha256_reaction_rules
        self.sha256_matching_rules = sha256_matching_rules

    def __repr__(self) -> str:
        return f"<Rules: {len(self._reaction_rules)} reaction rules, {len(self._matching_rules)} matching rules>"

    def get_reaction_rules(self, group_names: list[str] | None = None) -> list[ReactionRule]:
        """
        Get the reaction rules.

        :return: list of reaction rules
        """
        if group_names is not None:
            # Filter reaction rules by group name
            self._logger.debug(f"retrieving reaction rules for groups: {group_names}")

            return [r for r in self._reaction_rules if any(g in r.groups for g in group_names)]

        return self._reaction_rules

    def get_matching_rules(self) -> list[MatchingRule]:
        """
        Get the matching rules sorted by priori

        :return: sorted list of matching rules
        """
        matching_rules = self._matching_rules

        return matching_rules

    def check_for_duplicates(self) -> None:
        """
        Checks for duplicate items in matching rules.

        :raises ValueError: if duplicate items are found
        """
        logger = logging.getLogger(LOGGER_NAME)

        errors_seen = 0

        seen: dict[str, MatchingRule] = {}
        for rule in tqdm(self._matching_rules, desc="checking for duplicate matching rules"):
            curr_rid = rule.rid
            curr_mol = rule.mol

            try:
                inchikey = mol_to_inchikey(curr_mol)
            except Exception as e:
                raise RuntimeError(f"Failed to make InChIKey for rule {curr_rid}") from e

            if inchikey in seen:
                prev: MatchingRule = seen[inchikey]
                prev_rid = prev.rid

                logger.warning(
                    "Duplicate matching rule detected:\n"
                    f"- First: {prev_rid} ({mol_to_smiles(deepcopy(prev.mol), remove_tags=True)})\n"
                    f"- Second: {curr_rid} ({mol_to_smiles(deepcopy(curr_mol), remove_tags=True)})\n"
                    f"- InChIKey: {inchikey}"
                )
                errors_seen += 1

            seen[inchikey] = rule

        if errors_seen > 0:
            # Ask user input to contine or not; don't throw error
            logger.warning(f"Found {errors_seen} duplicate matching rules.")
            user_input = input("Do you want to continue? (y/n): ")
            if user_input.lower() != "y":
                raise ValueError("Duplicate matching rules found, exiting.")
            else:
                logger.warning("Continuing despite duplicate matching rules.")
        else:
            logger.info("No duplicate matching rules found.")


def get_path_default_reaction_rules() -> str:
    """
    Get the default path to the reaction rules JSON file.

    :return: path to the default reaction rules JSON file
    """
    return files(retromol.data).joinpath("default_reaction_rules.yml")


def get_path_default_matching_rules() -> str:
    """
    Get the default path to the matching rules JSON file.

    :return: path to the default matching rules JSON file
    """
    return files(retromol.data).joinpath("default_matching_rules.yml")


def get_path_default_wave_config() -> str:
    """
    Get the default path to the wave configuration JSON file.

    :return: path to the default wave configuration JSON file
    """
    return files(retromol.data).joinpath("default_wave_config.yml")


def load_rules_from_files(
    path_reaction_rules: str,
    path_matching_rules: str,
) -> Rules:
    """
    Load rules from JSON files.

    :param path_reaction_rules: path to the reaction rules JSON file
    :param path_matching_rules: path to the matching rules JSON file
    :return: Rules object containing the loaded rules
    """
    reaction_rules_src = open(path_reaction_rules).read()
    matching_rules_src = open(path_matching_rules).read()
    sha256_reaction_rules = sha256_hex(reaction_rules_src)
    sha256_matching_rules = sha256_hex(matching_rules_src)

    with open(path_reaction_rules) as fo:
        reaction_rules_data = yaml.safe_load(fo)
        reaction_rules = [ReactionRule.from_json_serializable_dict(i, r) for i, r in enumerate(reaction_rules_data)]

    with open(path_matching_rules) as fo:
        matching_rules_data = yaml.safe_load(fo)
        matching_rules = [MatchingRule.from_json_serializable_dict(i, r) for i, r in enumerate(matching_rules_data)]

    return Rules(reaction_rules, matching_rules, sha256_reaction_rules, sha256_matching_rules)
