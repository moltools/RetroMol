"""Shared retrieval machinery for figures 2, 3, and 4: fingerprint-only, alignment-only,
and combined (fingerprint prefilter + alignment re-rank) compound retrieval for a BGC
query, plus whole-molecule Tanimoto (Tc) scoring of a retrieved candidate against a
ground-truth compound.

"combined" mirrors gui/src/server/routes/discovery.py's run_discovery_query exactly:
fingerprint-cosine nearest neighbors first (cheap, DuckDB-vectorized, see
RetroMolDuckDB.closest), then re-rank that shortlist by normalized alignment score.
"alignment_only" has no such prefilter, so it's only tractable over a bounded
candidate pool (see benchmark/config.yaml's retrieval.alignment_only_candidate_pool_size).
"""

import random
import time
from dataclasses import dataclass, field
from typing import Any

import duckdb
import numpy as np
from rdkit import Chem

from common import DiscoveryContext, exact_label_group_accessions, normalize_for_alignment
from retromol.chem.fingerprint import calculate_tanimoto_similarity, mol_to_morgan_fingerprint
from retromol_alignment.pairwise import align
from retromol_alignment.ranking import reorder_target_chains
from retromol_database.duckdb import RetroMolDuckDB


@dataclass(frozen=True)
class CompoundRecord:
    entry_id: str
    smiles: str | None
    primary_sequence: list[str]


@dataclass
class RankedCandidate:
    entry_id: str
    score: float  # cosine similarity (fingerprint_only) or normalized alignment pct/100 (alignment_only, combined)


@dataclass(frozen=True)
class BgcQuery:
    entry_id: str
    accession: str
    fingerprint: np.ndarray
    primary_sequence: list[str]
    ground_truth_ids: list[str]


def group_accessions_from_config(
    db: RetroMolDuckDB, links: dict[str, dict], mibig_highlight_classes: dict
) -> dict[str, set[str]]:
    """{group_key: {accession, ...}} for every group in benchmark/config.yaml's
    mibig_highlight_classes, resolved via exact_label_group_accessions -- an
    accession only joins a group if its representative bgc entry's biosynthetic_class
    label set is EXACTLY that group's `biosynthetic_class_labels` (see that config
    key's comment): a hybrid BGC with an extra, unlisted label doesn't spuriously
    join any group it's only partially a match for."""
    group_label_sets = {
        key: frozenset(group["biosynthetic_class_labels"]) for key, group in mibig_highlight_classes.items()
    }
    return exact_label_group_accessions(db.con, links, group_label_sets)


def build_eval_bgcs(
    db: RetroMolDuckDB,
    links: dict[str, dict],
    group_accessions: dict[str, set[str]],
) -> dict[str, list[BgcQuery]]:
    """One BgcQuery per MIBiG accession that has both a bgc entry and >=1 associated
    compound entry (its "ground truth" set for retrieval evaluation), grouped into
    "all" plus each group in `group_accessions` (e.g. polyketide, nonribosomal-peptide
    accessions; see figure2_retrieval.py/figure3_cross_modal.py for how those sets
    are built)."""
    eligible = {
        accession: ids
        for accession, ids in links.items()
        if ids["bgc_entry_ids"] and ids["compound_entry_ids"]
    }

    all_bgc_ids = sorted({eid for ids in eligible.values() for eid in ids["bgc_entry_ids"]})
    entries_by_id = {e.id: e for e in db.get_entries(all_bgc_ids)}

    queries: dict[str, BgcQuery] = {}
    for accession, ids in eligible.items():
        # A BGC accession normally maps to exactly one region/entry; if antiSMASH
        # produced more than one region for the same accession, use the first
        # (stable, sorted) one as this accession's query. One query per accession
        # keeps every group's "n" comparable to a BGC/product count, not a raw region
        # count.
        bgc_id = ids["bgc_entry_ids"][0]
        entry = entries_by_id.get(bgc_id)
        if entry is None or not entry.primary_sequence:
            continue
        queries[accession] = BgcQuery(
            entry_id=bgc_id,
            accession=accession,
            fingerprint=np.array(entry.fingerprint, dtype=np.float32),
            primary_sequence=list(entry.primary_sequence),
            ground_truth_ids=list(ids["compound_entry_ids"]),
        )

    groups: dict[str, list[BgcQuery]] = {"all": list(queries.values())}
    for key, accessions in group_accessions.items():
        groups[key] = [q for acc, q in queries.items() if acc in accessions]

    return groups


def load_compound_universe(db_path: str) -> list[CompoundRecord]:
    con = duckdb.connect(db_path, read_only=True)
    try:
        rows = con.execute(
            "SELECT id, raw, primary_sequence FROM entries WHERE type = 'compound'"
        ).fetchall()
    finally:
        con.close()
    return [CompoundRecord(entry_id=r[0], smiles=r[1], primary_sequence=list(r[2])) for r in rows]


def sample_candidate_pool(
    universe: list[CompoundRecord],
    pool_size: int | None,
    seed: int,
    must_include_ids: set[str],
) -> list[CompoundRecord]:
    """Random (seeded) subsample of the compound universe, always keeping every entry
    in `must_include_ids` (e.g. this evaluated set's own ground-truth compounds), so
    recall@k stays measurable even when the pool is much smaller than the full DB."""
    if pool_size is None or pool_size >= len(universe):
        return universe

    must_have = [c for c in universe if c.entry_id in must_include_ids]
    rest = [c for c in universe if c.entry_id not in must_include_ids]
    n_random = max(0, pool_size - len(must_have))
    rng = random.Random(seed)
    sampled = rng.sample(rest, min(n_random, len(rest)))
    return must_have + sampled


def self_alignment_score(sequence: list[str], ctx: DiscoveryContext) -> float:
    score, _, _ = align(ctx.aligner, sequence, sequence, ctx.converter)
    return score


def normalized_pct(align_score: float, denom: float) -> float:
    if denom <= 0:
        return 0.0
    return max(0.0, min(1.0, align_score / denom))


def score_candidate_alignment(
    ctx: DiscoveryContext,
    query_seq_norm: list[str],
    self_score: float,
    score_mode: str,
    candidate: CompoundRecord,
) -> float | None:
    target = [normalize_for_alignment(n, ctx) for n in candidate.primary_sequence]
    if not target or not all(t in ctx.alphabet for t in target):
        return None

    _assignment_score, reordered_target, _inverted = reorder_target_chains(
        query_seq_norm, target, ctx.aligner, ctx.converter
    )
    flat_score, _, _ = align(ctx.aligner, query_seq_norm, reordered_target, ctx.converter)

    if score_mode == "longest_sequence":
        target_self = self_alignment_score(target, ctx)
        denom = max(self_score, target_self)
    else:
        denom = self_score

    return normalized_pct(flat_score, denom)


def alignment_only_rank(
    ctx: DiscoveryContext,
    query_seq_norm: list[str],
    self_score: float,
    candidates: list[CompoundRecord],
    score_mode: str,
    top_k: int,
) -> list[RankedCandidate]:
    scored: list[RankedCandidate] = []
    for cand in candidates:
        pct = score_candidate_alignment(ctx, query_seq_norm, self_score, score_mode, cand)
        if pct is not None:
            scored.append(RankedCandidate(cand.entry_id, pct))
    scored.sort(key=lambda r: r.score, reverse=True)
    return scored[:top_k]


def fingerprint_only_rank(db: RetroMolDuckDB, query_fp: np.ndarray, limit: int) -> list[RankedCandidate]:
    hits = db.closest(query_fp, limit=limit, entry_type="compound")
    return [RankedCandidate(h.entry.id, h.similarity) for h in hits]


def combined_rank(
    db: RetroMolDuckDB,
    ctx: DiscoveryContext,
    query_fp: np.ndarray,
    query_seq_norm: list[str],
    self_score: float,
    score_mode: str,
    prefilter_n: int,
    top_k: int,
) -> list[RankedCandidate]:
    shortlist = db.closest(query_fp, limit=prefilter_n, entry_type="compound")
    scored: list[RankedCandidate] = []
    for hit in shortlist:
        cand = CompoundRecord(hit.entry.id, hit.entry.raw, hit.entry.primary_sequence)
        pct = score_candidate_alignment(ctx, query_seq_norm, self_score, score_mode, cand)
        if pct is not None:
            scored.append(RankedCandidate(cand.entry_id, pct))
    scored.sort(key=lambda r: r.score, reverse=True)
    return scored[:top_k]


def timed_rank(fn, *args, **kwargs) -> tuple[list[RankedCandidate], float]:
    start = time.perf_counter()
    result = fn(*args, **kwargs)
    elapsed = time.perf_counter() - start
    return result, elapsed


class TanimotoCache:
    """Whole-molecule Morgan-fingerprint Tanimoto (Tc) between two compounds' SMILES,
    memoized per compound (each compound's own fingerprint is reused across every
    pair it appears in)."""

    def __init__(self, radius: int, num_bits: int, use_chirality: bool) -> None:
        self.radius = radius
        self.num_bits = num_bits
        self.use_chirality = use_chirality
        self._fp_cache: dict[str, Any] = {}

    def _fp(self, smiles: str | None):
        if smiles is None:
            return None
        if smiles not in self._fp_cache:
            mol = Chem.MolFromSmiles(smiles)
            self._fp_cache[smiles] = (
                mol_to_morgan_fingerprint(mol, self.radius, self.num_bits, self.use_chirality)
                if mol is not None
                else None
            )
        return self._fp_cache[smiles]

    def tc(self, smiles_a: str | None, smiles_b: str | None) -> float | None:
        fp_a, fp_b = self._fp(smiles_a), self._fp(smiles_b)
        if fp_a is None or fp_b is None:
            return None
        return calculate_tanimoto_similarity(fp_a, fp_b)
