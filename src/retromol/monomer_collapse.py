"""Collapse monomers into structural (and optionally name-based) groups, deterministically."""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field

from retromol.chem import (
    ExplicitBitVect,
    Mol,
    calc_tanimoto_similarity_rdkit,
    ecfp4,
    mol_to_inchikey,
    mol_to_smiles,
    standardize_from_smiles,
)
from retromol.helpers import blake64_hex


def inchikeys(mol: Mol) -> tuple[str, str]:
    """
    Get the InChIKeys for a molecule.

    :param mol: input molecule
    :return: tuple of (full InChIKey, connectivity InChIKey)
    """
    ik_full = mol_to_inchikey(mol)
    ik_conn = ik_full.split("-")[0]
    return ik_full, ik_conn


@dataclass
class Monomer:
    """
    Normalized view of each input record, kept for deterministic sorting and lookup.

    :param idx: original index in input list
    :param name: monomer name
    :param input_smiles: original input SMILES string
    :param mol: standardized RDKit Mol object (or None if standardization failed)
    :param can_smi: canonical SMILES string (or None if mol is None)
    :param ik_full: full InChIKey (or None if mol is None)
    :param ik_conn: connectivity InChIKey (or None if mol is None)
    :param fp: ECFP4 fingerprint (or None if mol is None)
    """

    idx: int
    name: str
    input_smiles: str
    mol: Mol | None = None
    can_smi: str | None = None
    ik_full: str | None = None
    ik_conn: str | None = None
    fp: ExplicitBitVect | None = None


@dataclass
class Group:
    """
    A group of monomers collapsed either by structure or by explicit name.

    :param gid: group ID
    :param rep_idx: index of the representative monomer
    :param members: list of member monomer indices
    :param token_fine: 64-bit hex over canonical SMILES (no stereo) OR name (for name-groups)
    :param rep_can_smi: canonical SMILES for the representative
    :param kind: "struct" or "name"
    :param name_key: the name string used for name-based groups
    """

    gid: int
    rep_idx: int  # index of the representative monomer
    members: list[int] = field(default_factory=list)
    kind: str = "struct"  # "struct" or "name"
    token_fine: str = ""  # 64-bit hex over canonical SMILES (struct) OR name (for name-groups)
    rep_can_smi: str = ""  # canonical SMILES for the representative
    name_key: str | None = None  # the name string used for name-based groups

    token_coarse: str = ""  # e.g., scaffold hash or family hash (placeholder for now)
    rep_scaffold_smi: str = ""  # if needed later for computing scaffolds, safe default


@dataclass
class NameSimilarityConfig:
    """
    Configure similarity bits among name-collapsed groups.

    :param family_of: maps a group name -> family string (None to skip family)
    :param family_weight: weight for each family token (per distinct member name)
    :param pair_weight: base weight multiplied by pairwise similarity in [0, 1]
    :param pairwise: sparse matrix (dict of dict) of explicit similarities. Example: {'serine': {'homoserine': 0.85}}
    :param symmetric: if True, treat pairwise as symmetric (use max(a->b, b->a))
    :param family_repeat_scale: integer scaling factor for family tokens
    :param pair_repeat_scale: integer scaling factor for pairwise tokens
    :param ancestors_of: maps a group name -> ancestor string (None to skip ancestor)
    :param ancestor_repeat_scale: integer scaling factor for ancestor tokens
    """

    family_of: Callable[[str], str | None] | None = None
    family_weight: float = 0.30
    pair_weight: float = 0.60
    pairwise: Mapping[str, Mapping[str, float]] = field(default_factory=dict)
    symmetric: bool = True
    family_repeat_scale: int = 2
    pair_repeat_scale: int = 2
    ancestors_of: Callable[[str], str | None] | None = None
    ancestor_repeat_scale: int = 0


class DSU:
    """Disjoint Set Union (Union-Find) data structure for efficient component merging."""

    def __init__(self, n: int) -> None:
        """
        Initialize DSU with n elements (0 to n-1).

        :param n: number of elements
        """
        self.p = list(range(n))
        self.r = [0] * n

    def find(self, x: int) -> int:
        """
        Find the representative of the set containing x, with path compression.

        :param x: element to find
        :return: representative element of the set
        """
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]  # path halving
            x = self.p[x]
        return x

    def union(self, a: int, b: int) -> None:
        """
        Union the sets containing elements a and b.

        :param a: first element
        :param b: second element
        """
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.r[ra] < self.r[rb]:
            self.p[ra] = rb
        elif self.r[ra] > self.r[rb]:
            self.p[rb] = ra
        else:
            self.p[rb] = ra
            self.r[ra] += 1


def collapse_monomers_order_invariant(
    records: Iterable[tuple[str, str]],
    keep_stereo: bool = False,
    tanimoto_thresh: float = 0.85,
    collapse_by_name: Iterable[str] | None = None,
) -> tuple[list[Group], list[Monomer]]:
    """
    Deterministic grouping independent of input order (but still RDKit/version dependent).

    Algorithm:
    1) Normalize each (name, SMILES) -> Monomer (std Mol, can_smi, InChIKeys, ECFP4)
    2) Split indices into name-driven vs structure-driven pools
    3) Structural pool:
        a) Exact union by full-IK and by connectivity-IK
        b) Block by rough size (bitcount bucker) and union pairs with Tanimoto >= threshold
    4) Emit groups deterministically:
        a) Name groups in sorted (name) with members sorted by stable key
        b) Structural components by representative stable key
    5) Return monomers sorted by stable key for reproducibility

    :param records: iterable of (name, SMILES) tuples for monomers
    :param keep_stereo: whether to retain stereochemistry during standardization
    :param tanimoto_thresh: Tanimoto similarity threshold for structural grouping
    :param collapse_by_name: optional iterable of names to always collapse by name
    :return: tuple of (list of Groups, list of Monomers)
    """
    collapse_set = set(collapse_by_name or [])

    # Build Monomer table (skip invalid SMILES unless in name-collapsed)
    monomers: list[Monomer] = []
    for i, (name, smi) in enumerate(records):
        mol = standardize_from_smiles(smi, keep_stereo=keep_stereo) if smi else None
        if mol is None and name not in collapse_set:
            continue
        can_smi = mol_to_smiles(mol, isomeric=keep_stereo, canonical=True) if mol is not None else None
        ik_full = ik_conn = None
        fp = None
        if mol is not None:
            ik_full, ik_conn = inchikeys(mol)
            fp = ecfp4(mol)
        monomers.append(
            Monomer(
                idx=i,
                name=name,
                input_smiles=smi or "",
                mol=mol,
                can_smi=can_smi,
                ik_full=ik_full,
                ik_conn=ik_conn,
                fp=fp,
            )
        )

    # Stable key used globally to kill order effects
    def mkey(m: Monomer) -> tuple[str, str, str, int]:
        """
        Stable key for monomer sorting and representative selection.

        :param m: Monomer object
        :return: tuple key
        """
        return (m.can_smi or "", m.name or "", m.input_smiles or "", m.idx)

    # Helper: get Monomer by original idx quickly
    by_idx: dict[int, Monomer] = {m.idx: m for m in monomers}

    # Partition into name vs structural pools
    name_idxs = [m.idx for m in monomers if m.name in collapse_set]
    struct_idxs = [m.idx for m in monomers if m.name not in collapse_set]

    # Dense indexing for DSU only over structural pool
    struct_pos = {idx: pos for pos, idx in enumerate(sorted(struct_idxs))}
    dsu = DSU(len(struct_pos))

    # Tier 1: exact unions by InChIKeys (full, then connectivity)
    by_full: dict[str, list[int]] = defaultdict(list)
    by_conn: dict[str, list[int]] = defaultdict(list)
    for i in struct_idxs:
        m = by_idx[i]
        if m.ik_full:
            by_full[m.ik_full].append(i)
        if m.ik_conn:
            by_conn[m.ik_conn].append(i)

    for bucket in list(by_full.values()) + list(by_conn.values()):
        # Deterministic chaining unions after sorting
        if len(bucket) >= 2:
            sb = sorted(bucket)
            for a, b in zip(sb[:-1], sb[1:], strict=True):
                dsu.union(struct_pos[a], struct_pos[b])

    # Tier 2: similarity unions inside coarse bitcount blocks
    def bitcount(fp: ExplicitBitVect | None) -> int:
        """
        Count the number of bits set in an RDKit ExplicitBitVect fingerprint.

        :param fp: RDKit ExplicitBitVect fingerprint
        :return: Number of bits set
        """
        return int(fp.GetNumOnBits()) if fp is not None else 0

    blocks: dict[tuple[str, int], list[int]] = defaultdict(list)
    for i in struct_idxs:
        m: Monomer = by_idx[i]
        # Single "no scaffold" channel: bucket by rough size (per 16 bits set)
        blocks[("<NOSCAF>", bitcount(m.fp) // 16)].append(i)

    # Deterministic within each block: sort by mkey, then scan upper triangle
    for _, bucket in sorted(blocks.items(), key=lambda kv: kv[0]):
        if len(bucket) < 2:
            continue
        bucket = sorted(bucket, key=lambda i: mkey(by_idx[i]))
        for ai in range(len(bucket)):
            ma = by_idx[bucket[ai]]
            if ma.fp is None:
                continue
            for bi in range(ai + 1, len(bucket)):
                mb = by_idx[bucket[bi]]
                if mb.fp is None:
                    continue
                if calc_tanimoto_similarity_rdkit(ma.fp, mb.fp) >= tanimoto_thresh:
                    dsu.union(struct_pos[ma.idx], struct_pos[mb.idx])

    # Emit groups deterministically
    groups: list[Group] = []

    # Name groups: emit in sorted (name) order; members sorted by mkey
    names_sorted = sorted({by_idx[i].name for i in name_idxs})
    for nm in names_sorted:
        mems = [i for i in name_idxs if by_idx[i].name == nm]
        mems_sorted = sorted(mems, key=lambda i: mkey(by_idx[i]))
        rep = by_idx[mems_sorted[0]]
        groups.append(
            Group(
                gid=len(groups),
                rep_idx=rep.idx,
                members=mems_sorted,
                token_fine=blake64_hex(f"NAME:{nm}"),
                rep_can_smi=rep.can_smi or "",
                kind="name",
                name_key=nm,
            )
        )

    # Structural components: gather, choose representative by mkey, and sort components by rep key
    comps: dict[int, list[int]] = defaultdict(list)
    for i in struct_idxs:
        root = dsu.find(struct_pos[i])
        comps[root].append(i)

    comp_infos: list[tuple[tuple[str, str, str, int], list[int]]] = []
    for comp in comps.values():
        comp_sorted = sorted(comp, key=lambda i: mkey(by_idx[i]))
        rep = by_idx[comp_sorted[0]]
        rep_key = (rep.can_smi or "", rep.name or "", rep.idx)
        comp_infos.append((rep_key, comp_sorted))

    for _, comp_sorted in sorted(comp_infos, key=lambda t: t[0]):
        rep = by_idx[comp_sorted[0]]
        groups.append(
            Group(
                gid=len(groups),
                rep_idx=rep.idx,
                members=comp_sorted,
                token_fine=blake64_hex(rep.can_smi or ""),
                rep_can_smi=rep.can_smi or "",
                kind="struct",
            )
        )

    # Return monomers in a determinisitic order as well
    monomers_sorted = sorted(monomers, key=mkey)
    return groups, monomers_sorted


def tokens_for_groups(
    groups: list[Group],
    weight_fine: float = 1.0,
    weight_coarse: float = 0.4,
    extra_bags: Iterable[dict[str, float]] | None = None,
) -> dict[str, float]:
    """
    Build a determinisitc token bag (fine/coarse), sorting to keep stable accumulation order.

    :param groups: list of Group objects
    :param weight_fine: weight for fine tokens (name or canonical SMILES)
    :param weight_coarse: weight for coarse tokens (scaffold-level, not implemented here)
    :param extra_bags: optional coarse tokens (e.g., scaffold/family) if present in Group
    :return: dictionary mapping tokens to their accumulated weights
    """
    bag: dict[str, float] = {}

    # Stable accumulation order
    for g in sorted(groups, key=lambda g: (g.kind, g.rep_can_smi, g.gid)):
        if g.token_fine:
            prefix = "N:" if g.kind == "name" else "F:"
            key = f"{prefix}{g.token_fine}"
            bag[key] = bag.get(key, 0.0) + weight_fine

        if g.token_coarse:
            keyc = f"C:{g.token_coarse}"
            bag[keyc] = bag.get(keyc, 0.0) + weight_coarse

    if extra_bags:
        for eb in extra_bags:
            for k, v in eb.items():
                bag[k] = bag.get(k, 0.0) + float(v)

    return bag


def assign_to_existing_groups(
    smi: str,
    groups: list[Group],
    monomers: list[Monomer],
    keep_stereo: bool = False,
    tanimoto_thresh: float = 0.85,
) -> int | None:
    """
    Deterministically assign SMILES into existing structural groups.

    :param smi: input SMILES string.
    :param groups: list of existing Group objects.
    :param monomers: list of Monomer objects corresponding to the groups.
    :param keep_stereo: whether to retain stereochemistry during standardization.
    :param tanimoto_thresh: Tanimoto similarity threshold for structural grouping.
    :return: group ID if assigned, else None.
    """
    mol = standardize_from_smiles(smi, keep_stereo=keep_stereo)
    if mol is None:
        return None

    ik_full, ik_conn = inchikeys(mol)
    fp_new = ecfp4(mol)

    # Representatives, determinisitically sorted
    reps = [(g.gid, monomers[g.rep_idx]) for g in groups if g.kind == "struct"]
    reps.sort(key=lambda t: (t[1].can_smi or "", t[0]))

    fullIK_to_gid = {m.ik_full: gid for gid, m in reps if m.ik_full}
    connIK_to_gid = {m.ik_conn: gid for gid, m in reps if m.ik_conn}

    # Exact InChIKey
    if ik_full in fullIK_to_gid:
        return fullIK_to_gid[ik_full]

    # Connectivity InChIKey
    if ik_conn in connIK_to_gid:
        return connIK_to_gid[ik_conn]

    # Tanimoto similarity fallback
    best_gid, best_sim = None, 0.0
    for gid, rep in reps:
        if rep.fp is None:
            continue
        sim = calc_tanimoto_similarity_rdkit(fp_new, rep.fp)
        if sim > best_sim:
            best_gid, best_sim = gid, sim

    return best_gid if best_sim >= tanimoto_thresh else None
