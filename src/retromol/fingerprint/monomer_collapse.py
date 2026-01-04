"""Collapse monomers into structural (and optionally name-based) groups, deterministically."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable

from rdkit.Chem.rdchem import Mol
from rdkit.DataStructs.cDataStructs import ExplicitBitVect

from retromol.chem.mol import mol_to_smiles, standardize_from_smiles
from retromol.chem.fingerprint import mol_to_morgan_fingerprint, calculate_tanimoto_similarity
from retromol.model.rules import MatchingRule
from retromol.utils.hashing import blake64_hex


@dataclass(frozen=True, slots=True)
class Monomer:
    """
    Normalized monomer view derived from a MatchingRule. Only keep what we need for determinisic grouping.

    :var rid: rule ID (index in matching rules)
    :var name: monomer name
    :var mol: RDKit Mol object
    :var can_smi: canonical SMILES (isomeric if specified)
    :var fp: ECFP4 fingerprint
    """

    rid: int
    name: str
    mol: Mol
    can_smi: str
    fp: ExplicitBitVect

    @classmethod
    def from_matching_rule(
        cls,
        matching_rule: MatchingRule,
        keep_stereo: bool = False,
        morgan_radius: int = 2,
        morgan_num_bits: int = 2048,
    ) -> "Monomer":
        """
        Create Monomer from MatchingRule.

        :param matching_rule: MatchingRule object
        :return: Monomer instance
        """
        canonical_smiles = mol_to_smiles(matching_rule.mol, isomeric=keep_stereo, canonical=True)
        morgan_fingerprint = mol_to_morgan_fingerprint(matching_rule.mol, radius=morgan_radius, num_bits=morgan_num_bits)

        return cls(
            rid=matching_rule.id,
            name=matching_rule.name,
            mol=matching_rule.mol,
            can_smi=canonical_smiles,
            fp=morgan_fingerprint,
        )


@dataclass
class Group:
    """
    Deterministic structural group.

    :var gid: group ID
    :var rep_rid: representative Monomer.rid
    :var members: tuple of member Monomer.rid values (sorted)
    :var token: stable token based on rep_can_smi or name
    :var rep_can_smi: representative canonical SMILES
    """

    gid: int
    rep_rid: int
    members: tuple[int, ...]
    token: str
    rep_can_smi: str


def collapse_monomers(
    matching_rules: Iterable[MatchingRule],
    keep_stereo: bool = False,
    tanimoto_threshold: float = 0.6,
    morgan_radius: int = 2,
    morgan_num_bits: int = 2048,
) -> tuple[list[Group], list[Monomer]]:
    """
    Deterministic grouping independent of input order (but still RDKit/version dependent).

    :param records: iterable of (name, SMILES) tuples for monomers
    :param keep_stereo: whether to retain stereochemistry during standardization
    :param tanimoto_threshold: Tanimoto similarity threshold for structural grouping
    :param morgan_radius: radius for Morgan fingerprint
    :param morgan_num_bits: number of bits for Morgan fingerprint
    :return: tuple of (list of Groups, list of Monomers)
    """
    # Build Monomer table from matching rules
    monomers: list[Monomer] = [
        Monomer.from_matching_rule(
            rl,
            keep_stereo=keep_stereo,
            morgan_radius=morgan_radius,
            morgan_num_bits=morgan_num_bits,
        )
        for rl in matching_rules
    ]

    # Stable order independent of input order
    monomers.sort(key=lambda m: (m.can_smi, m.rid))

    # Excat groups by canonical smiles
    # key: can_smi -> list of monomer rids
    by_smi: dict[str, list[Monomer]] = {}
    for m in monomers:
        by_smi.setdefault(m.can_smi, []).append(m)

    # Sort exact groups by (can_smi) so deterministic
    exact_reps: list[Monomer] = []
    exact_members: list[list[Monomer]] = []
    for can_smi in sorted(by_smi.keys()):
        ms = sorted(by_smi[can_smi], key=lambda x: x.rid)
        exact_reps.append(ms[0])  # deterministic representative
        exact_members.append(ms)

    # Similarity collapse across exact groups
    # We collapse exact groups into "structural families" deterministically
    rep_indices: list[int] = []  # indices into exact_reps that are final representatives
    assigned_to: list[int] = [-1] * len(exact_reps)  # exact-group i -> rep index (index in rep_indices list)

    for i, rep_i in enumerate(exact_reps):
        # Find first earlier representative that matches by similarity
        found_rep_slot = None
        for slot, rep_group_idx in enumerate(rep_indices):
            rep_j = exact_reps[rep_group_idx]
            sim = calculate_tanimoto_similarity(rep_i.fp, rep_j.fp)
            if sim >= tanimoto_threshold:
                found_rep_slot = slot
                break

        if found_rep_slot is None:
            # New representative
            rep_indices.append(i)
            assigned_to[i] = len(rep_indices) - 1
        else:
            # Assign to existing representative
            assigned_to[i] = found_rep_slot

    # Emit final groups deterministically
    groups: list[Group] = []
    for gid, rep_group_idx in enumerate(rep_indices):
        rep = exact_reps[rep_group_idx]

        # Gather members from all exact-groups assigned to this rep slot
        member_rids: list[int] = []
        for i in range(len(exact_reps)):
            if assigned_to[i] != gid:
                continue
            member_rids.extend(m.rid for m in exact_members[i])

        member_rids = sorted(set(member_rids))  # unique and sorted
        groups.append(Group(
            gid=gid,
            rep_rid=rep.rid,
            members=tuple(member_rids),
            token=blake64_hex(rep.can_smi),
            rep_can_smi=rep.can_smi,
        ))

    return groups, monomers


def assign_to_group(
    smiles: str,
    groups: list[Group],
    monomers: list[Monomer],
    keep_stereo: bool = False,
    tanimoto_threshold: float = 0.6,
    morgan_radius: int = 2,
    morgan_num_bits: int = 2048,
) -> Group | None:
    """
    Assign a new monomer (by SMILES) to an existing group if similar enough.

    :param smiles: SMILES string of new monomer 
    :param groups: existing groups
    :param monomers: existing monomers
    :param keep_stereo: whether to retain stereochemistry during standardization
    :param tanimoto_threshold: Tanimoto similarity threshold for structural grouping
    :param morgan_radius: radius for Morgan fingerprint
    :param morgan_num_bits: number of bits for Morgan fingerprint
    :return: assigned Group or None if no match
    """
    mol = standardize_from_smiles(smiles, keep_stereo=keep_stereo)
    if mol is None:
        raise ValueError(f"could not standardize SMILES {smiles} for group assignment")
    
    # Create gid -> group mapping
    gid_to_group = {g.gid: g for g in groups}

    # Get canonical SMILES and fingerprint for SMILES-to-assign
    can_smi = mol_to_smiles(mol, include_tags=False, canonical=True, isomeric=keep_stereo)
    fp_new = mol_to_morgan_fingerprint(mol, radius=morgan_radius, num_bits=morgan_num_bits)

    # Exact canonical SMILES -> group (fast path)
    rep_smi_to_gid = {g.rep_can_smi: g.gid for g in groups}
    gid = rep_smi_to_gid.get(can_smi)
    if gid is not None:
        return gid_to_group[gid]
    
    # Similarity fallback vs. representative fingerprints
    # Build rid -> Monomer lookup once (monomers contain fp + can_smi)
    by_rid = {m.rid: m for m in monomers}

    best_gid: int | None = None
    best_sim: float = -1.0

    # Deterministic iteration: sort by gid
    for g in sorted(groups, key=lambda x: x.gid):
        rep = by_rid.get(g.rep_rid)
        if rep is None:
            continue  # should not happen if monomers list matches groups

        sim = calculate_tanimoto_similarity(fp_new, rep.fp)
        if sim > best_sim:
            best_sim = sim
            best_gid = g.gid

    if best_gid is None:
        return None
    
    if best_sim >= tanimoto_threshold:
        return gid_to_group[best_gid]
    
    return None
