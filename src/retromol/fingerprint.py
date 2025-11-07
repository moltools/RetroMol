"""This module contains functions for generating hashed fingerprints from k-mers."""

import hashlib
import os
import pickle
import re
import struct
from collections.abc import Callable, Iterable, Sequence
from datetime import datetime
from itertools import combinations
from typing import Any, TypeVar

import numpy as np
import yaml
from numpy.typing import NDArray

from retromol.graph import iter_kmers
from retromol.helpers import blake64_hex, sha256_hex
from retromol.io import Result
from retromol.monomer_collapse import (
    Group,
    NameSimilarityConfig,
    assign_to_existing_groups,
    collapse_monomers_order_invariant,
)
from retromol.readout import mapping_to_graph, optimal_mappings
from retromol.rules import MatchingRule

T = TypeVar("T")


def _norm_token(tok: object, none_sentinel: str = "<NONE>") -> bytes:
    """
    Turn a token (possibly None) into stable bytes.

    :param tok: token to normalize (str, int, float, or None)
    :param none_sentinel: string to use for None tokens
    :return: bytes representation of the token
    """
    if tok is None:
        return none_sentinel.encode("utf-8")

    # Strings/ints are common; fall back to repr for others
    if isinstance(tok, (str, int, float)):
        return str(tok).encode("utf-8")

    return repr(tok).encode("utf-8")


def _family_token(fam: str) -> str:
    """
    Generate a family token.

    :param fam: family name
    :return: family token string
    """
    fam = fam or ""
    return f"NF:{blake64_hex(f'FAM:{fam.lower()}')}"


def _pair_token(a: str, b: str) -> str:
    """
    Generate a pairwise token for two names.

    :param a: first name
    :param b: second name
    :return: pairwise token string
    """
    a, b = sorted([(a or "").lower(), (b or "").lower()])
    return f"NS:{blake64_hex(f'PAIR:{a}|{b}')}"


def _hash_kmer_tokens(
    tokens: Sequence[bytes],
    n_bits: int,
    n_hashes: int,
    seed: int = 0,
    k_salt: int = 0,
) -> list[int]:
    """
    Map a tokenized k-mer (as bytes) to n_hashes bit indices in [0, n_bits).

    :param tokens: sequence of bytes tokens (e.g. from _norm_token)
    :param n_bits: number of bits in the fingerprint
    :param n_hashes: number of hash indices to produce
    :param seed: global seed for hashing
    :param k_salt: salt value specific to the k-mer length (to decorrelate lengths)
    :return: list of bit indices

    .. note:: Deterministic across runs/machines. Different k values get a salt.
    """
    data = b"\x1f".join(tokens)  # unit separator

    idxs: list[int] = []
    for i in range(n_hashes):
        # Include both global seed and per-hash index, plus a per-k salt
        salted = data + struct.pack(">III", seed, i, k_salt)
        digest = hashlib.blake2b(salted, digest_size=8).digest()
        val = int.from_bytes(digest, "big") % n_bits
        idxs.append(val)

    return idxs


def kmers_to_fingerprint(
    kmers: Iterable[Sequence[Any]],
    n_bits: int = 2048,
    n_hashes_per_kmer: int | Callable[[int], int] = 2,
    seed: int = 42,
    none_policy: str = "keep",
    counted: bool = False,
    count_dtype: Any = np.uint32,
) -> NDArray[np.generic]:
    """
    Build a hashed fingerprint from an iterable of tokenized k-mers.

    :param kmers: iterable of k-mers, where each k-mer is a sequence of tokens (str, int, float, or None)
    :param n_bits: number of bits in the fingerprint
    :param n_hashes_per_kmer: number of hash indices to produce per k-mer (int or callable that takes k-mer length
        as input and returns the number of hashes).
    :param seed: global seed for hashing.
    :param none_policy: policy for handling None tokens: "keep" (treat as a special token), "skip-token"
        (omit the token), or "drop-kmer" (skip the entire k-mer).
    :param counted: if True, produce a count vector instead of a binary vector
    :param count_dtype: data type for counts (if counted is True)
    :return: fingerprint as a numpy array of shape (n_bits,)
    """
    if n_bits <= 0:
        raise ValueError("n_bits must be positive")

    # Normalize n_hashes_per_kmer to callable
    if isinstance(n_hashes_per_kmer, int):
        if n_hashes_per_kmer <= 0:
            raise ValueError("n_hashes_per_kmer must be positive")

        def _nh(_: int) -> int:
            return n_hashes_per_kmer
    else:
        _nh: Callable[[int], int] = n_hashes_per_kmer

    # Allocate output
    if counted:
        vec = np.zeros(n_bits, dtype=count_dtype)
    else:
        vec = np.zeros(n_bits, dtype=np.uint8)

    # Main loop
    for kmer in kmers:
        if none_policy == "drop-kmer" and any(t is None for t in kmer):
            continue

        # Normalize per token
        normd: list[bytes] = []
        for t in kmer:
            if t is None:
                if none_policy == "skip-token":
                    continue
                normd.append(_norm_token(None))
            else:
                normd.append(_norm_token(t))
        if not normd:
            continue

        n_hashes = _nh(len(kmer))
        if n_hashes <= 0:
            continue

        # Simple salt tied to (normalized) k-mer length
        k_salt = len(normd)

        idxs = _hash_kmer_tokens(
            normd,
            n_bits=n_bits,
            n_hashes=n_hashes,
            seed=seed,
            k_salt=k_salt,
        )

        if counted:
            # Increment counts; duplicates in idxs will accumulate
            vec[idxs] += 1
        else:
            # Set bits to 1 (binary)
            vec[idxs] = 1

    return vec


def cosine_similarity(fp1: NDArray[np.int8], fp2: NDArray[np.int8]) -> float:
    """
    Cosine similarity for fingerprints.

    :param fp1: first fingerprint (1D array)
    :param fp2: second fingerprint (1D array)
    :return: cosine similarity in [0, 1]
    """
    a = np.asarray(fp1)
    b = np.asarray(fp2)

    # Ensure 1D
    a = a.ravel()
    b = b.ravel()
    if a.shape != b.shape:
        raise ValueError(f"Different lengths: {a.shape} vs {b.shape}")

    # Upcast to float to avoid integer overflow and match sklearn
    a = a.astype(np.float64, copy=False)
    b = b.astype(np.float64, copy=False)

    # Compute cosine
    dot = float(np.dot(a, b))
    na = float(np.linalg.norm(a))
    nb = float(np.linalg.norm(b))

    if na == 0.0 or nb == 0.0:
        return 0.0

    return dot / (na * nb)


def tanimoto_similarity(fp1: NDArray[np.int8], fp2: NDArray[np.int8]) -> float:
    """
    Tanimoto similarity for molecular fingerprints (binary or count-based).

    :param fp1: first fingerprint (1D array)
    :param fp2: second fingerprint (1D array)
    :return: Tanimoto similarity in [0, 1]
    """
    a = np.asarray(fp1)
    b = np.asarray(fp2)

    # Ensure 1D
    a = a.ravel()
    b = b.ravel()
    if a.shape != b.shape:
        raise ValueError(f"Different lengths: {a.shape} vs {b.shape}")

    # Upcast to float to prevent overflow and ensure precision
    a = a.astype(np.float64, copy=False)
    b = b.astype(np.float64, copy=False)

    # Dot product = intersection term
    ab = float(np.dot(a, b))
    aa = float(np.dot(a, a))
    bb = float(np.dot(b, b))

    denom = aa + bb - ab
    if denom == 0.0:
        return 0.0
    return ab / denom


def get_kmers(seq: tuple[T, ...], k: int) -> list[tuple[T, ...]]:
    """
    Return all contiguous, bidirectional k-mers (subtuples of length k) from a tuple.

    :param seq: input sequence as a tuple
    :param k: k-mer length
    :return: list of k-mers (as tuples)
    """
    if k <= 0:
        return []
    n = len(seq)
    if k > n:
        return []
    forward_kmers = [seq[i : i + k] for i in range(n - k + 1)]
    backward_kmers = [tuple(reversed(kmer)) for kmer in forward_kmers]
    return forward_kmers + backward_kmers


class FingerprintGenerator:
    """Class to generate fingerprints based on monomer collapse groups."""

    def __init__(
        self,
        matching_rules_yaml: str | None,
        keep_stereo: bool = False,
        tanimoto_threshold: float = 0.85,
        collapse_by_name: list[str] | None = None,
        name_similarity: NameSimilarityConfig | None = None,
    ) -> None:
        """
        Initialize FingerprintGenerator.

        :param matching_rules_yaml: path to matching rules YAML file
        :param keep_stereo: whether to retain stereochemistry during standardization
        :param tanimoto_threshold: Tanimoto similarity threshold for structural grouping
        :param collapse_by_name: optional list of names to always collapse by name
        :param name_similarity: optional configuration for name similarity
        :raises FileNotFoundError: if the matching rules YAML file does not exist
        """
        if not os.path.exists(matching_rules_yaml):
            raise FileNotFoundError(f"Matching rules YAML file not found: {matching_rules_yaml}")

        # Load matching rules and compute SHA256 hash for provenance
        with open(matching_rules_yaml) as f:
            matching_rules_src = f.read()
            sha256_matching_rules = sha256_hex(matching_rules_src)
            self.sha256_matching_rules = sha256_matching_rules

            # Parse out the matching rules, and turn into records
            matching_rules_data = yaml.safe_load(matching_rules_src)
            records = []
            for i, rule_data in enumerate(matching_rules_data):
                matching_rule = MatchingRule.from_json_serializable_dict(i, rule_data)
                records.append((matching_rule.rid, matching_rule.smiles))

        # Collapse monomers into groups
        groups, monomers = collapse_monomers_order_invariant(
            records,
            keep_stereo=keep_stereo,
            tanimoto_thresh=tanimoto_threshold,
            collapse_by_name=collapse_by_name,
        )
        self.collapse_by_name = collapse_by_name or []
        self.groups = groups
        self.monomers = monomers
        self.name_similarity = name_similarity

        self.keep_stereo = keep_stereo
        self.tanimoto_threshold = tanimoto_threshold

        # For speedup
        self._assign_cache: dict[tuple[str | None, str], Group | None] = {}
        self._token_bytes_cache: dict[object, bytes] = {}

    def __repr__(self) -> str:
        """
        String representation of the FingerprintGenerator.

        :return: string representation
        """
        return (
            f"FingerprintGenerator(num_groups={len(self.groups)}, "
            f"num_monomers={len(self.monomers)}, "
            f"keep_stereo={self.keep_stereo}, "
            f"tanimoto_threshold={self.tanimoto_threshold}, "
            f"collapse_by_name={self.collapse_by_name}, "
            f"name_similarity={self.name_similarity})"
        )

    def assign_to_group(self, smi: str, name: str | None = None) -> Group:
        """
        Assign a new monomer to an existing group based on its SMILES.

        :param name: name of the monomer
        :param smi: SMILES string of the monomer
        :return: group ID if assigned, None otherwise
        """
        # Cache key: only use name when we're collapsing by that name
        key = (name if (name is not None and name in self.collapse_by_name) else None, smi)

        # Return from cache (including cached None) if present
        g = self._assign_cache.get(key)
        if g is not None or key in self._assign_cache:
            return g  # may be None

        # Name branch: cache the hit when found
        if key[0] is not None:
            # Deterministically scan existing roups
            for gg in self.groups:
                if gg.kind == "name" and gg.name_key == key[0]:
                    self._assign_cache[key] = gg
                    return gg
            # If we intended to collapse by name but no such group exists, that's an error
            raise ValueError(f"No existing name-based group found for name: {name}")

        # Structure branch: assign based on Tanimoto similarity
        gid = assign_to_existing_groups(
            smi=smi,
            groups=self.groups,
            monomers=self.monomers,
            keep_stereo=self.keep_stereo,
            tanimoto_thresh=self.tanimoto_threshold,
        )

        # gid can be 0; only None means "not assigned"
        g = self.groups[gid] if gid is not None else None

        # Cache result (including None) so we don't recompute on repeats
        self._assign_cache[key] = g
        return g

    def fingerprint_from_result(
        self,
        result: Result,
        num_bits: int = 2048,
        kmer_sizes: list[int] | None = None,
        kmer_weights: dict[int, int] | None = None,
        strict: bool = True,
        counted: bool = False,
    ) -> NDArray[np.int8] | None:
        """
        Generate a fingerprint from a RetroMolResult.

        :param result: RetroMol Result object
        :param num_bits: number of bits in the fingerprint
        :param kmer_sizes: list of k-mer sizes to consider
        :param kmer_weights: weights for each k-mer size. Determines how many bits each k-mer sets.
        :param strict: if True, verify that the matching rules SHA256 matches.
        :param counted: if True, count the number of times each k-mer appears.
        :return: fingerprint as a numpy array, or None if no monomers found.
        """
        if strict:
            sha256_self = self.sha256_matching_rules
            sha256_result = result.sha256_matching_rules
            if sha256_self != sha256_result:
                raise ValueError(
                    "Mismatch in matching rules SHA256: FingerprintGenerator was "
                    "created with different matching rules than the Result"
                )

        # Defaults
        if kmer_sizes is None:
            kmer_sizes = [1, 2, 3]
        if kmer_weights is None:
            kmer_weights = {1: 16, 2: 4, 3: 2}

        # Resolve similarity config
        cfg = self.name_similarity
        family_of = cfg.family_of if cfg and cfg.family_of is not None else (lambda n: None)
        pairwise = (cfg.pairwise if cfg else {}) or {}
        symmetric = bool(cfg.symmetric) if cfg else True
        fam_rep = max(0, int(cfg.family_repeat_scale)) if cfg else 0
        pair_rep = max(0, int(cfg.pair_repeat_scale)) if cfg else 0
        ancestors_of = getattr(cfg, "ancestors_of", None) if cfg else None
        anc_rep = max(0, int(getattr(cfg, "ancestor_repeat_scale", 0))) if cfg else 0

        # Gather optimal mappings
        oms = [om for om in optimal_mappings(result)]

        # Get tagged SMILES from result
        tagged_smiles = result.get_input_smiles(remove_tags=False)

        fps = []
        for om in oms:
            om_graph = mapping_to_graph(tagged_smiles, om)

            token_kmers: list[tuple[str, ...]] = []
            names_per_kmer: list[list[str]] = []
            sizes_per_kmer: list[int] = []

            for kmer_size in kmer_sizes:
                for kmer in iter_kmers(om_graph, kmer_size):
                    tokenized_kmer = []
                    names_in_kmer = []

                    for node in kmer:
                        node_data = om_graph.nodes[node]
                        node_id = node_data.get("identity")
                        node_smiles = node_data.get("smiles")
                        if node_smiles is None:
                            raise ValueError("Node in mapping graph missing 'smiles' attribute")

                        group = self.assign_to_group(node_smiles, name=node_id)
                        token = group.token_fine if group is not None else None
                        tokenized_kmer.append(token)

                        if node_id is not None:
                            names_in_kmer.append(node_id)

                    token_kmers.append(tuple(tokenized_kmer))
                    names_per_kmer.append(names_in_kmer)
                    sizes_per_kmer.append(kmer_size)

            # Inject similarity "virtual 1-mers" (families and pairwise), per k-mer
            if cfg and (fam_rep > 0 or pair_rep > 0):
                # Family tokens: once per name in the k-mer, repeated fam_rep times
                if fam_rep > 0:
                    for names in names_per_kmer:
                        for nm in sorted(set(n for n in names if n)):
                            fam_val = family_of(nm)
                            # Accept str, iterable of str, or None
                            # This allows for multiple families per name
                            if fam_val is None:
                                families = []
                            elif isinstance(fam_val, (list, tuple, set)):
                                families = list(fam_val)
                            else:
                                families = [fam_val]

                            for fam in sorted({f for f in families if f}):
                                ftok = _family_token(fam)
                                if not ftok:
                                    continue
                                for _ in range(fam_rep):
                                    token_kmers.append((ftok,))

                # Pairwise tokens: for unordered name pairs in the same k-mer
                if pair_rep > 0 and pairwise:
                    for names in names_per_kmer:
                        uniq = sorted(set(n for n in names if n))
                        if len(uniq) < 2:
                            continue
                        for a, b in combinations(uniq, 2):
                            s = float(pairwise.get(a, {}).get(b, 0.0))
                            if symmetric:
                                s = max(s, float(pairwise.get(b, {}).get(a, 0.0)))
                            if s <= 0.0:
                                continue
                            reps = int(round(s * pair_rep))
                            if reps <= 0:
                                continue
                            ptoken = _pair_token(a, b)
                            for _ in range(reps):
                                token_kmers.append((ptoken,))

            # Ancestor supertokens
            if ancestors_of and anc_rep > 0:
                # Small local helper: stable ancestor token w/ level namespace
                def _anc_tok(level: int, anc: str) -> str:
                    anc = (anc or "").lower()
                    return f"AN:{level}:{blake64_hex(f'ANC:{level}:{anc}')}"

                # Ancestor 1-mers: for each name, emit all ancestors in its path
                for names in names_per_kmer:
                    for nm in set(n for n in names if n):
                        path = ancestors_of(nm) or []  # e.g., ["polyketide", "polyketide_type_A", "A1"]
                        for lvl, anc in enumerate(path):
                            tok = _anc_tok(lvl, anc)
                            for _ in range(anc_rep):
                                token_kmers.append((tok,))

                # Ancestor k-mers: for each window, for each ancestor level present at all positions
                for names, ksize in zip(names_per_kmer, sizes_per_kmer, strict=True):
                    if ksize <= 1:
                        continue
                    pos_paths = [(ancestors_of(nm) or []) if nm else [] for nm in names]
                    if not pos_paths or any(len(p) == 0 for p in pos_paths):
                        # Require every position to have at least one ancestor (root level)
                        continue
                    max_depth = min(len(p) for p in pos_paths)  # only levels common to all positions
                    for lvl in range(max_depth):
                        # Form one ancestor k-mer at this level by taking the ancestor token per position
                        kmer_tok = tuple(_anc_tok(lvl, pos_paths[i][lvl]) for i in range(ksize))
                        for _ in range(anc_rep):
                            token_kmers.append(kmer_tok)

            # Hash default + virtual kmers
            fp = kmers_to_fingerprint(
                token_kmers,
                n_bits=num_bits,
                n_hashes_per_kmer=lambda k: kmer_weights.get(k, 1),
                seed=42,
                none_policy="skip-token",
                counted=counted,
            )
            fps.append(fp)

        # Stack fingerprints from all optimal mappings as rows
        if not fps:
            return None

        stacked_fps = np.stack(fps, axis=0)  # shape (num_mappings, num_bits)
        return stacked_fps

    def fingerprint_from_kmers(
        self,
        kmers: Iterable[Sequence[tuple[str | None, str | None]]],
        num_bits: int = 2048,
        kmer_weights: dict[int, int] | None = None,
        *,
        counted: bool = False,
        none_policy: str = "skip-token",
        allow_raw_name_token: bool = True,
        raise_on_unknown_named_group: bool = True,
    ) -> NDArray[np.int8]:
        """
        Build a fingerprint directly from user-provided k-mers of (name, smiles) tuples.

        Each item in a k-mer must be a 2-tuple: (name|None, smiles|None).
        At least one element of the tuple must be non-None.

        Name-only tokens (e.g. ('chlorination', None)):
          - If `name` is present and in `self.collapse_by_name`, we resolve to the
            existing name-based group token (same as used elsewhere).
            If such a group doesn't exist and `raise_on_unknown_named_group` is True,
            a ValueError is raised; otherwise we fall back to a stable name token
            (when `allow_raw_name_token` is True).
          - If `name` is present but NOT in `collapse_by_name`, we optionally
            use a stable name token if `allow_raw_name_token` is True.

        Structure-only tokens (e.g. (None, 'CCO')):
          - Assigned structurally via `assign_to_group(smiles, name=None)`;
            if no structural group is found, the item contributes no token.

        Mixed tokens (e.g. ('A12', 'CCO')):
          - If `name` is in `collapse_by_name`, we resolve via the name group.
            Otherwise we assign structurally; `name` is still passed to
            `assign_to_group` (it can help for caches, but does not force name grouping).

        :param kmers: iterable of k-mers; each k-mer is a sequence of (name, smiles) tuples
        :param num_bits: size of the fingerprint
        :param kmer_weights: hash multiplicity per k (default {1:16, 2:4, 3:2})
        :param counted: if True, produce a count vector; otherwise binary
        :param none_policy: how to handle None tokens at the token level
            ('keep', 'skip-token', 'drop-kmer'). Defaults to 'skip-token'
        :param allow_raw_name_token: if True, when a name isn't in `collapse_by_name`
            (or no name-group exists), use a stable name token so it still hashes
        :param raise_on_unknown_named_group: if True and a name is in `collapse_by_name`
            but no name-based group exists, raise ValueError
        :return: 1D numpy array (num_bits,) with the fingerprint
        """
        if kmer_weights is None:
            kmer_weights = {1: 16, 2: 4, 3: 2}

        def _stable_name_token(nm: str) -> str:
            # Deterministic, case-insensitive token for raw names
            return f"NM:{blake64_hex('NAME:' + (nm or '').lower())}"

        token_kmers: list[tuple[str | None, ...]] = []
        names_per_kmer: list[list[str]] = []
        sizes_per_kmer: list[int] = []

        for kmer in kmers:
            if not kmer:
                continue

            toks: list[str | None] = []
            names_here: list[str] = []

            for item in kmer:
                if not isinstance(item, tuple) or len(item) != 2:
                    raise TypeError("Each k-mer item must be a (name|None, smiles|None) tuple")
                name, smi = item
                if name is None and smi is None:
                    # Explicitly ignore completely empty items
                    continue

                tok: str | None = None

                # Priority 1: explicit collapse-by-name
                if name is not None and name in (self.collapse_by_name or []):
                    try:
                        # assign_to_group will find the existing name-group without needing a real SMILES
                        g = self.assign_to_group(smi="", name=name)
                    except ValueError:
                        # No such name group exists
                        if raise_on_unknown_named_group:
                            raise
                        g = None
                    if g is not None:
                        tok = g.token_fine
                    elif allow_raw_name_token and name is not None:
                        tok = _stable_name_token(name)

                # Priority 2: structure-based assignment
                elif smi:
                    g = self.assign_to_group(smi=smi, name=name)
                    tok = g.token_fine if g is not None else None

                # Priority 3: raw name token
                elif allow_raw_name_token and name is not None:
                    tok = _stable_name_token(name)

                # If tok remains None, the item contributes nothing
                toks.append(tok)

                if name:
                    names_here.append(name)

            # Emit this k-mer (even if some items were None; none_policy will handle)
            if toks:
                token_kmers.append(tuple(toks))
                names_per_kmer.append(names_here)
                sizes_per_kmer.append(len(toks))  # logical k size

        # Inject similarity virtual tokens
        cfg = self.name_similarity
        if cfg:
            family_of = cfg.family_of if cfg.family_of is not None else (lambda n: None)
            pairwise = cfg.pairwise or {}
            symmetric = bool(cfg.symmetric) if cfg.symmetric is not None else True
            fam_rep = max(0, int(cfg.family_repeat_scale or 0))
            pair_rep = max(0, int(cfg.pair_repeat_scale or 0))
            ancestors_of = getattr(cfg, "ancestors_of", None)
            anc_rep = max(0, int(getattr(cfg, "ancestor_repeat_scale", 0)))

            # helper for stable ancestor token
            def _anc_tok(level: int, anc: str) -> str:
                anc = (anc or "").lower()
                return f"AN:{level}:{blake64_hex(f'ANC:{level}:{anc}')}"

            # Families
            if fam_rep > 0:
                for names in names_per_kmer:
                    for nm in sorted(set(n for n in names if n)):
                        fam_val = family_of(nm)
                        families = (
                            []
                            if fam_val is None
                            else (list(fam_val) if isinstance(fam_val, (list, tuple, set)) else [fam_val])
                        )
                        for fam in sorted({f for f in families if f}):
                            ftok = _family_token(fam)
                            for _ in range(fam_rep):
                                token_kmers.append((ftok,))

            # Pairwise
            if pair_rep > 0 and pairwise:
                for names in names_per_kmer:
                    uniq = sorted(set(n for n in names if n))
                    if len(uniq) < 2:
                        continue
                    for a, b in combinations(uniq, 2):
                        s = float(pairwise.get(a, {}).get(b, 0.0))
                        if symmetric:
                            s = max(s, float(pairwise.get(b, {}).get(a, 0.0)))
                        reps = int(round(max(0.0, s) * pair_rep))
                        if reps > 0:
                            ptoken = _pair_token(a, b)
                            for _ in range(reps):
                                token_kmers.append((ptoken,))

            # Ancestors (1-mers and aligned k-mers)
            if ancestors_of and anc_rep > 0:
                # 1-mers
                for names in names_per_kmer:
                    for nm in set(n for n in names if n):
                        path = ancestors_of(nm) or []
                        for lvl, anc in enumerate(path):
                            tok = _anc_tok(lvl, anc)
                            for _ in range(anc_rep):
                                token_kmers.append((tok,))
                # aligned k-mers by ancestor level
                for names, ksize in zip(names_per_kmer, sizes_per_kmer, strict=True):
                    if ksize <= 1:
                        continue
                    pos_paths = [(ancestors_of(nm) or []) if nm else [] for nm in names]
                    if not pos_paths or any(len(p) == 0 for p in pos_paths):
                        continue
                    max_depth = min(len(p) for p in pos_paths)
                    for lvl in range(max_depth):
                        kmer_tok = tuple(_anc_tok(lvl, pos_paths[i][lvl]) for i in range(ksize))
                        for _ in range(anc_rep):
                            token_kmers.append(kmer_tok)

        if not token_kmers:
            # Return an all-zero vector of the requested type/shape to keep behavior predictable
            return kmers_to_fingerprint([], n_bits=num_bits, n_hashes_per_kmer=1, counted=counted)

        fp = kmers_to_fingerprint(
            token_kmers,
            n_bits=num_bits,
            n_hashes_per_kmer=lambda k: kmer_weights.get(k, 1),
            seed=42,
            none_policy=none_policy,
            counted=counted,
        )
        return fp

    def save(self, path: str) -> None:
        """
        Serialize this FingerprintGenerator to a binary file.

        :param path: output file path
        .. note:: uses pickle; not secure against untrusted sources
        """
        payload = {
            "__format__": "retromol.FingerprintGenerator",
            "__version__": 1,
            "__created__": datetime.utcnow().isoformat() + "Z",
            # Core config/state
            "keep_stereo": self.keep_stereo,
            "tanimoto_threshold": self.tanimoto_threshold,
            "collapse_by_name": self.collapse_by_name,
            "name_similarity": self.name_similarity,
            "sha256_matching_rules": self.sha256_matching_rules,
            # Precomputed data
            "groups": self.groups,
            "monomers": self.monomers,
        }
        with open(path, "wb") as f:
            pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, path: str) -> "FingerprintGenerator":
        """
        Load a FingerprintGenerator previously saved with .save().

        :param path: input file path
        :return: FingerprintGenerator instance
        .. note:: uses pickle; not secure against untrusted sources
        """
        with open(path, "rb") as f:
            payload = pickle.load(f)

        # Basic validation
        fmt = payload.get("__format__")
        ver = int(payload.get("__version__", -1))
        if fmt != "retromol.FingerprintGenerator" or ver != 1:
            raise ValueError(f"Unrecognized FingerprintGenerator save format/version: {fmt} v{ver}")

        # Build instance without invoking __init__
        self = cls.__new__(cls)

        # Restore core config/state
        self.keep_stereo = bool(payload["keep_stereo"])
        self.tanimoto_threshold = float(payload["tanimoto_threshold"])
        self.collapse_by_name = list(payload["collapse_by_name"]) if payload["collapse_by_name"] else []
        self.name_similarity = payload["name_similarity"]
        self.sha256_matching_rules = payload["sha256_matching_rules"]

        # Restore precomputed data
        self.groups = payload["groups"]
        self.monomers = payload["monomers"]

        # Recreate transient caches
        self._assign_cache = {}
        self._token_bytes_cache = {}

        return self


def polyketide_family_of(name: str) -> list[str] | None:
    """
    Simple polyketide family extractor based on name pattern.

    :param name: monomer name
    :return: [family, subfamily] or None if not a polyketide
    """
    n = (name or "").strip()
    if not n:
        return None
    is_polyketide = re.match(r"^[ABCD]\d+$", n) is not None
    if is_polyketide:
        family = "polyketide"
        subfamily = n[0]
        return [family, subfamily]
    return None


def polyketide_ancestors_of(name: str) -> list[str]:
    """
    Simple polyketide ancestor extractor based on name pattern.

    :param name: monomer name
    :return: list of ancestors (e.g., ["polyketide", "polyketide_type_A", "A1"])
    """
    n = (name or "").strip().upper()
    if re.match(r"^[ABCD]$", n):
        return ["polyketide", f"polyketide_type_{n}"]
    if re.match(r"^[ABCD]\d+$", n):
        return ["polyketide", f"polyketide_type_{n[0]}", n]
    return []
