"""This module contains functions for generating hashed fingerprints from k-mers."""

import hashlib
import struct
import logging
from typing import Any, Literal, Callable, Iterable, Sequence, Iterator, TypeVar

import numpy as np
from numpy.typing import NDArray

_BIOCRACKER = False
try:
    from biocracker.query.modules import (
        LinearReadout as BioCrackerLinearReadout,
        Module,
        NRPSModule,
        PKSModule,
        PKSExtenderUnit,
    )
    _BIOCRACKER = True
except ImportError:
    pass

from retromol.model.result import Result
from retromol.model.rules import MatchingRule
from retromol.model.reaction_graph import MolNode
from retromol.utils.hashing import blake64_hex

from retromol.fingerprint.monomer_collapse import Group, collapse_monomers, assign_to_group


log = logging.getLogger(__name__)


NonePolicy = Literal["keep", "skip-token", "drop-kmer"]


_MISS = object()
DEFAULT_KMER_WEIGHTS: dict[int, int] = {1: 1, 2: 1}
DEFAULT_KMER_SIZES: list[int] = [1, 2]


T = TypeVar("T")
Direction = Literal["forward", "backward", "both"]


def iter_kmers_sequence(items: Sequence[T], k: int, direction: Direction = "both") -> Iterator[tuple[T, ...]]:
    """
    Iterate over k-mers of size k from a sequence of items.

    :param items: sequence of items
    :param k: size of the k-mers
    :param direction: direction of k-mers to generate: "forward", "backward", or "both"
    :return: iterator over k-mers as tuples
    """
    if k < 1:
        raise ValueError("k must be at least 1")
    
    n = len(items)
    if k > n:
        return  # no kmers possible
    
    # Special case: 1-mers have no meaningful direction
    if k == 1:
        for x in items:
            yield (x,)
        return
    
    def forward() -> Iterator[tuple[T, ...]]:
        """
        Yield forward k-mers.
        """
        for i in range(0, n - k + 1):
            yield tuple(items[i : i + k])
    
    def backward() -> Iterator[tuple[T, ...]]:
        """
        Yield backward k-mers.
        """
        for km in forward():
            yield km[::-1]

    if direction == "forward":
        yield from forward()
    elif direction == "backward":
        yield from backward()
    elif direction == "both":
        yield from forward()
        yield from backward()
    else:
        raise ValueError(f"invalid direction: {direction}")


def encode_family_token(fam: str) -> str:
    """
    Generate a family token.

    :param fam: family name
    :return: family token string
    """
    fam = fam or ""
    return f"NF:{blake64_hex(f'FAM:{fam.lower()}')}"


def normalize_token(tok: object, none_sentinel: str = "<NONE>") -> bytes:
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


def hash_kmer_tokens(
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
    num_bits: int = 2048,
    num_hashes_per_kmer: int | Callable[[int], int] = 2,
    seed: int = 42,
    none_policy: NonePolicy = "keep",
    counted: bool = False,
) -> NDArray[np.generic]:
    """
    Build a hashed fingerprint from an iterable of tokenized k-mers.

    :param kmers: iterable of k-mers, where each k-mer is a sequence of tokens (str, int, float, or None)
    :param num_bits: number of bits in the fingerprint
    :param num_hashes_per_kmer: number of hash indices to produce per k-mer (int or callable that takes k-mer length
        as input and returns the number of hashes).
    :param seed: global seed for hashing.
    :param none_policy: policy for handling None tokens: "keep" (treat as a special token), "skip-token"
        (omit the token), or "drop-kmer" (skip the entire k-mer).
    :param counted: if True, produce a count vector instead of a binary vector
    :return: fingerprint as a numpy array of shape (n_bits,)
    """
    if num_bits <= 0:
        raise ValueError("n_bits must be positive")

    # Normalize n_hashes_per_kmer to callable
    if isinstance(num_hashes_per_kmer, int):
        if num_hashes_per_kmer <= 0:
            raise ValueError("num_hashes_per_kmer must be positive")

        def _nh(_: int) -> int:
            return num_hashes_per_kmer
    else:
        _nh: Callable[[int], int] = num_hashes_per_kmer

    # Allocate output
    if counted:
        fp = np.zeros(num_bits, dtype=np.uint32)
    else:
        fp = np.zeros(num_bits, dtype=np.uint8)

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
                normd.append(normalize_token(None))
            else:
                normd.append(normalize_token(t))
        if not normd:
            continue

        n_hashes = _nh(len(kmer))
        if n_hashes <= 0:
            continue

        # Simple salt tied to (normalized) k-mer length
        k_salt = len(normd)

        idxs = hash_kmer_tokens(
            normd,
            n_bits=num_bits,
            n_hashes=n_hashes,
            seed=seed,
            k_salt=k_salt,
        )

        if counted:
            # Increment counts; duplicates in idxs will accumulate
            fp[idxs] += 1
        else:
            # Set bits to 1 (binary)
            fp[idxs] = 1

    return fp 


class FingerprintGenerator:
    """
    Class to generate fingerprints based on monomer collapse groups.
    """

    def __init__(
        self,
        matching_rules: Iterable[MatchingRule],
        keep_stereo: bool = False,
        tanimoto_threshold: float = 0.6,
        morgan_radius: int = 2,
        morgan_num_bits: int = 2048,
    ) -> None:
        """
        Initialize FingerprintGenerator.

        :param matching_rules: iterable of MatchingRule objects for monomer identification
        :param keep_stereo: whether to keep stereochemistry when collapsing monomers
        :param tanimoto_threshold: Tanimoto similarity threshold for collapsing monomers
        :param morgan_radius: radius for Morgan fingerprinting when collapsing monomers
        :param morgan_num_bits: number of bits for Morgan fingerprinting when collapsing monomers
        """
        matching_rules = list(matching_rules)

        groups, monomers = collapse_monomers(
            matching_rules,
            keep_stereo=keep_stereo,
            tanimoto_threshold=tanimoto_threshold,
            morgan_radius=morgan_radius,
            morgan_num_bits=morgan_num_bits,
        )

        self.groups = groups
        self.monomers = monomers

        self.keep_stereo = keep_stereo
        self.tanimoto_threshold = tanimoto_threshold
        self.morgan_radius = morgan_radius
        self.morgan_num_bits = morgan_num_bits

        # For speedup
        self._assign_cache: dict[str, Group | None] = {}
        self._token_bytes_cache: dict[object, bytes] = {}

    def assign_to_group(self, smiles: str) -> Group | None:
        """
        Assign a new monomer to an existing group based on its SMILES.

        :param smiles: SMILES string of the monomer
        :return: assigned Group or None if no match
        """
        # SMILES was checked before; return from cache
        cached = self._assign_cache.get(smiles, _MISS)
        if cached is not _MISS:
            return cached  # can be group or None

        # Structure branch: assign based on Tanimoto similarity
        group = assign_to_group(
            smiles=smiles,
            groups=self.groups,
            monomers=self.monomers,
            keep_stereo=self.keep_stereo,
            tanimoto_threshold=self.tanimoto_threshold,
            morgan_radius=self.morgan_radius,
            morgan_num_bits=self.morgan_num_bits,
        )

        # Cache result (including None) so we don't recompute on repeats
        self._assign_cache[smiles] = group
    
        return group
    
    def ancestor_list_for_node(self, node: MolNode) -> list[str | None]:
        """
        Return full ancestor hierarchy for a node.

        :param node: MolNode to get ancestors for
        :return: list of ancestor tokens (str or None)
        """
        anc: list[str] = []

        if node.is_identified and node.identity.matched_rule.ancestor_tokens:
            anc.extend(node.identity.matched_rule.ancestor_tokens)

        return anc

    def fingerprint_from_result(
        self,
        result: Result,
        num_bits: int = 2048,
        kmer_sizes: list[int] | None = None,
        kmer_weights: dict[int, int] | None = None,
        counted: bool = False,
    ) -> NDArray[np.int8]:
        """
        Generate a fingerprint from a RetroMolResult.

        :param result: RetroMol Result object
        :param num_bits: number of bits in the fingerprint
        :param kmer_sizes: list of k-mer sizes to consider
        :param kmer_weights: weights for each k-mer size. Determines how many bits each k-mer sets.
        :param counted: if True, count the number of times each k-mer appears.
        :return: fingerprint as a numpy array
        """
        # Default kmer_sizes
        if kmer_sizes is None:
            kmer_sizes = DEFAULT_KMER_SIZES

        # Default kmer_weights
        if kmer_weights is None:
            kmer_weights = DEFAULT_KMER_WEIGHTS

        # Retrieve AssemblyGraph from Result
        a = result.linear_readout.assembly_graph

        # Calculate kmers from AssemblyGraph
        tokenized_kmers: list[tuple[str | None, ...]] = []

        for kmer_size in kmer_sizes:
            for kmer in a.iter_kmers(k=kmer_size):

                # Ancestral tokens for items in kmer
                per_item_ancestors: list[list[str | None]] = []
                
                for item in kmer:
                    # Structural token is the lowest level ancestor
                    ancestors: list[str | None] = []

                    # First get the structural token (lowest level ancestor)
                    if item.is_identified:
                        smiles = item.smiles
                        ancestors.append(g.token if (g := self.assign_to_group(smiles)) is not None else None)
                    else:
                        ancestors.append(None)

                    # Then get the rest of the ancestors
                    # We reverse to have the highest level ancestor last
                    ancestors.extend(reversed(self.ancestor_list_for_node(item)))

                    per_item_ancestors.append(ancestors)
                
                assert len(per_item_ancestors) == len(kmer), "length mismatch in ancestor tokens"

                # Get tokenized kmer from every level of ancestor
                max_depth = max(len(anc) for anc in per_item_ancestors)
                for level in range(max_depth):
                    tokenized_kmers.append(tuple(
                        anc[level] if level < len(anc) else None
                        for anc in per_item_ancestors
                    ))

        # Gather additional 1-mer virtual family tokens (defined in matching rules); only once per found monomer
        for node in a.monomer_nodes():
            ident = node.identity if node.is_identified else None
            if ident is not None:
                for fam_tok in ident.matched_rule.family_tokens:
                    tokenized_kmers.append((encode_family_token(fam_tok),))

        # Hash kmers
        fp = kmers_to_fingerprint(
            tokenized_kmers,
            num_bits=num_bits,
            num_hashes_per_kmer=lambda k: kmer_weights.get(k, 1),
            seed=42,
            none_policy="keep",
            counted=counted,
        )

        return fp
    
    def fingerprint_from_biocracker_readout(
        self,
        readout: BioCrackerLinearReadout,
        by_orf: bool = False,
        num_bits: int = 2048,
        kmer_sizes: list[int] | None = None,
        kmer_weights: dict[int, int] | None = None,
        counted: bool = False,
    ) -> NDArray[np.int8]:
        """
        Generate a fingerprint from a BioCracker LinearReadout.

        :param readout: BioCracker LinearReadout object
        :param by_orf: if True, generate fingerprint per ORF instead of per region
        :param num_bits: number of bits in the fingerprint
        :param kmer_sizes: list of k-mer sizes to consider
        :param kmer_weights: weights for each k-mer size. Determines how many bits each k-mer sets.
        :param counted: if True, count the number of times each k-mer appears.
        :return: fingerprint as a numpy array
        :raises ImportError: if biocracker is not installed
        :raises ValueError: if unsupported module type is encountered
        :raises AssertionError: if length mismatches occur in token lists
        """
        if not _BIOCRACKER:
            raise ImportError("biocracker is not installed; cannot generate fingerprint from biocracker readout")

        # Default kmer_sizes
        if kmer_sizes is None:
            kmer_sizes = DEFAULT_KMER_SIZES

        # Default kmer_weights
        if kmer_weights is None:
            kmer_weights = DEFAULT_KMER_WEIGHTS

        # Calculate kmers from BioCracker's linear readout
        tokenized_kmers: list[tuple[str | None, ...]] = []

        ordered = readout.biosynthetic_order(by_orf=by_orf)

        # Normalize to an interable of (readout_id, modules_in_readout)
        if by_orf:
            # ordered: list[tuple[str, list[Module]]]
            readout_blocks: list[tuple[str, list[Module]]] = list(ordered)
        else:
            # ordered: list[Module] -> treat as a single block
            readout_blocks = [(readout.id, list(ordered))]
        
        for _, modules in readout_blocks:
            for kmer_size in kmer_sizes:
                for kmer in iter_kmers_sequence(modules, k=kmer_size, direction="both"):
                
                    # Ancestral tokens for items in kmer
                    per_module_ancestors: list[list[str | None]] = []

                    for module in kmer:
                        # Structural token is the lowest level ancestor
                        ancestors: list[str | None] = []

                        if isinstance(module, NRPSModule):
                            # Extract SMILES of predicted substrate
                            if module.substrate is not None:
                                if module.substrate.name == "graminine":
                                    # Graminine SMILES was incorrect in BioCracker versions <2.0.1; fix here for backwards compatibility
                                    smiles = r"O=NN(O)CCC[C@H](N)(C(=O)O)"
                                else:
                                    smiles = module.substrate.smiles
                            else:
                                smiles = None

                            if smiles is not None:
                                ancestors.append(g.token if (g := self.assign_to_group(smiles)) is not None else None)
                            else:
                                # No predicted substrate
                                ancestors.append(None) 

                            # We don't add ancestral tokens for NRPSModule

                        elif isinstance(module, PKSModule):
                            # PKSModule has no structural token
                            ancestors.append(None)

                            # Extract ancestral tokens
                            match module.substrate.extender_unit:
                                case PKSExtenderUnit.PKS_A: ancestors.extend(["A", "PKS"])
                                case PKSExtenderUnit.PKS_B: ancestors.extend(["B", "PKS"])
                                case PKSExtenderUnit.PKS_C: ancestors.extend(["C", "PKS"])
                                case PKSExtenderUnit.PKS_D: ancestors.extend(["D", "PKS"])
                        
                        else:
                            # Unsupported module type
                            log.warning(f"Unsupported module type: {type(module)}")
                            ancestors.append(None)

                        per_module_ancestors.append(ancestors)

                    assert len(per_module_ancestors) == len(kmer), "length mismatch in ancestor tokens"
                    
                    # Get tokenized kmer from every level of ancestor
                    max_depth = max(len(anc) for anc in per_module_ancestors)
                    for level in range(max_depth):
                        tokenized_kmers.append(tuple(
                            anc[level] if level < len(anc) else None
                            for anc in per_module_ancestors
                        ))

        # Add modifiers as family tokens
        for modifier in readout.modifiers:
            tokenized_kmers.append((encode_family_token(modifier),))

        # Hash kmers
        fp = kmers_to_fingerprint(
            tokenized_kmers,
            num_bits=num_bits,
            num_hashes_per_kmer=lambda k: kmer_weights.get(k, 1),
            seed=42,
            none_policy="keep",
            counted=counted,
        )

        return fp 
    