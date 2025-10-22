"""Module for matching reaction graph nodes to motifs."""

import logging
from collections.abc import Iterable
from copy import deepcopy
from typing import Any

from retromol import chem, config, rules


def match_mol_greedily(
    mol: chem.Mol, rls: list[rules.MatchingRule], sch: bool = False
) -> tuple[str, dict[str, Any]] | None:
    """
    Match a molecule to a motif.

    :param mol: RDKit molecule to match
    :param rls: list of matching rules (motifs)
    :param sch: whether to match stereochemistry
    :return: tuple of motif ID and properties if matched, else None
    .. note:: this function uses a greedy approach to match a molecule to a motif
    """
    for rl in rls:
        if rid := rl.is_match(mol, sch):
            return rid, rl.props

    return None


def greedy_max_set_cover(enc_to_mol: dict[int, chem.Mol], nodes: list[int]) -> list[int]:
    """
    Find biggest non-overlapping set of mol nodes in the reaction graph.

    :param enc_to_mol: mapping of encoding to RDKit molecule
    :param nodes: list of node encodings to consider for set cover
    :return: list of selected node encodings
    """
    # Create subsets of atom mappings per node.
    subsets: list[tuple[int, set[int]]] = list()
    for node in nodes:
        mol = enc_to_mol[node]
        tags = {atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetIsotope() != 0}
        subsets.append((node, tags))

    # Sort subsets by size of atom mappings, from largest to smallest.
    sorted_subsets = sorted(subsets, key=lambda x: len(x[1]), reverse=True)

    # Perform greedy set cover algorithm.
    selected_subsets: list[int] = []
    covered_elements: set[int] = set()
    for node, subset in sorted_subsets:
        uncovered_elements = subset - covered_elements

        # Make sure that only a subset is selected if all elements are uncovered.
        if uncovered_elements != subset:
            continue

        if uncovered_elements:
            selected_subsets.append(node)
            covered_elements.update(uncovered_elements)

    return selected_subsets


def solve_exact_cover_with_priority(
    enc_to_mol: dict[str, "chem.Mol"],
    nodes_A: list[str],
    nodes_B: list[str],
    required_tags: Iterable[int],
) -> tuple[list[int], list[int]]:
    """
    Partition `required_tags` into disjoint node-tag sets drawn from nodes_A ∪ nodes_B.

    Objective (lexicographic):
      1) Maximize the number of required tags covered by A (identified coverage)
      2) Among those, minimize the number of A nodes (prefer single big A over multiple small A)
      3) Then minimize total nodes

    :param enc_to_mol: mapping of encoding to RDKit molecule
    :param nodes_A: list of node encodings in set A (identified)
    :param nodes_B: list of node encodings in set B (unidentified)
    :param required_tags: iterable of required tags to cover
    :return: tuple of selected nodes from A and B
    :raises ValueError: if no exact cover exists
    """
    logger = logging.getLogger(config.LOGGER_NAME)

    req: set[int] = set(required_tags)

    def node_tags(node: str) -> set[int]:
        mol = enc_to_mol[node]
        return set(chem.get_tags_mol(mol))

    # Build candidate subsets (filter out nodes that contain tags outside the required set)
    candA = [(n, node_tags(n)) for n in nodes_A]
    candB = [(n, node_tags(n)) for n in nodes_B]
    candA = [(n, ts) for (n, ts) in candA if ts and ts.issubset(req)]
    candB = [(n, ts) for (n, ts) in candB if ts and ts.issubset(req)]

    # Quick impossibility check: every required tag must appear in at least one candidate
    tag_to_candidates: dict[int, list[tuple[str, str, set[int]]]] = {t: [] for t in req}
    for src, pool in (("A", candA), ("B", candB)):
        for n, ts in pool:
            for t in ts:
                tag_to_candidates[t].append((src, n, ts))

    for t in req:
        if not tag_to_candidates[t]:
            error_msg = f"No node covers required tag {t}; exact cover impossible."
            logger.error(error_msg)

            logger.error(f"Tags to cover: {req}")

            for node in nodes_A:
                mol = enc_to_mol[node]
                smiles_no_tags = chem.mol_to_smiles(deepcopy(mol), remove_tags=True)
                logger.debug(f"Node A {node}: SMILES {smiles_no_tags}")

            for node in nodes_B:
                mol = enc_to_mol[node]
                smiles_no_tags = chem.mol_to_smiles(deepcopy(mol), remove_tags=True)
                logger.debug(f"Node B {node}: SMILES {smiles_no_tags}")

            raise ValueError(error_msg)

    # Order candidates within each tag: A-first, larger sets first (helps reduce branching)
    for t in req:
        tag_to_candidates[t].sort(
            key=lambda x: (x[0] != "A", -len(x[2]))  # A before B; then bigger sets
        )

    # Greedy optimistic bound: how many tags can A at most still claim disjointly?
    def optimistic_A_tag_gain(remaining_tags: set[int], used_tags: set[int]) -> int:
        compat: list[tuple[str, set[int]]] = []
        for n, ts in candA:
            if ts and ts <= remaining_tags and ts.isdisjoint(used_tags):
                compat.append((n, ts))
        compat.sort(key=lambda x: -len(x[1]))  # larger tag-sets first
        covered: set[int] = set()
        for _, ts in compat:
            if ts.isdisjoint(covered):
                covered |= ts
                if covered >= remaining_tags:
                    break
        return len(covered)

    # Choose the "most constrained" remaining tag (fewest compatible candidates)
    def pick_most_constrained_tag(remaining: set[int], used: set[int]) -> int:
        best_t = None
        best_count = None
        for t in remaining:
            # Count only candidates compatible with used tags
            opts = 0
            for _, _, ts in tag_to_candidates[t]:
                if ts.isdisjoint(used) and ts <= remaining:
                    opts += 1
            if best_count is None or opts < best_count:
                best_t, best_count = t, opts
                if best_count == 1:
                    break  # can't get more constrained than 1
        return best_t

    best_solution: tuple[list[str], list[str]] | None = None
    best_A_tags = -1  # maximize
    best_A_nodes = 10**9  # minimize
    best_total = 10**9  # minimize

    def recurse(
        used_tags: set[int],
        used_A_tags: set[int],
        chosen_A: list[str],
        chosen_B: list[str],
    ):
        nonlocal best_solution, best_A_tags, best_A_nodes, best_total

        remaining = req - used_tags
        if not remaining:
            # Complete exact cover: compare by (|A-tags| desc, |A-nodes| asc, total asc)
            A_tags_cnt = len(used_A_tags)
            A_nodes_cnt = len(chosen_A)
            total = A_nodes_cnt + len(chosen_B)
            if (
                (A_tags_cnt > best_A_tags)
                or (A_tags_cnt == best_A_tags and A_nodes_cnt < best_A_nodes)
                or (A_tags_cnt == best_A_tags and A_nodes_cnt == best_A_nodes and total < best_total)
            ):
                best_solution = (chosen_A.copy(), chosen_B.copy())
                best_A_tags = A_tags_cnt
                best_A_nodes = A_nodes_cnt
                best_total = total
            return

        # Branch-and-bound: can we still beat the best A-tag count?
        # (Current |A-tags| + optimistic future A-tags) must exceed current best_A_tags.
        cur_A_tags = len(used_A_tags)
        if cur_A_tags + optimistic_A_tag_gain(remaining, used_tags) < best_A_tags:
            return

        # Pick the hardest tag next
        t = pick_most_constrained_tag(remaining, used_tags)

        # Try candidates that cover t, A-first (already ordered in tag_to_candidates)
        for src, n, ts in tag_to_candidates[t]:
            if not ts.isdisjoint(used_tags):
                continue  # not compatible (overlap)
            if not ts <= remaining:
                continue

            # Choose
            used_tags.update(ts)
            if src == "A":
                chosen_A.append(n)
                added_A = ts
                used_A_tags.update(added_A)
            else:
                chosen_B.append(n)
                added_A = None

            recurse(used_tags, used_A_tags, chosen_A, chosen_B)

            # Backtrack
            if src == "A":
                if added_A is None:
                    raise RuntimeError("Internal error: added_A should not be None for A nodes.")
                chosen_A.pop()
                used_A_tags.difference_update(added_A)
            else:
                chosen_B.pop()
            used_tags.difference_update(ts)

    recurse(set(), set(), [], [])

    if best_solution is None:
        raise ValueError("No exact cover found with the given A/B nodes for the required tags.")

    selA, selB = best_solution
    return selA, selB


def identify_nodes(
    encoding_to_mol: dict[str, chem.Mol],
    matching_rules: list[rules.MatchingRule],
    match_stereochemistry: bool = False,
) -> dict[str, dict[str, Any]]:
    """
    Identify nodes in a reaction graph that match given motifs.

    :param encoding_to_mol: mapping of encoding to RDKit molecule
    :param matching_rules: list of matching rules to apply
    :param match_stereochemistry: whether to match stereochemistry during identification
    :return: a dictionary mapping node encodings to motif IDs
    """
    # Try to identify nodes, keep those that match
    identity_mapping: dict[str, dict[str, Any]] = {}
    for node in encoding_to_mol.keys():
        mol = encoding_to_mol[node]
        matched = match_mol_greedily(mol, matching_rules, match_stereochemistry)
        if matched:
            rid, props = matched
            identity_mapping[node] = {"identity": rid, "props": props}

    return identity_mapping
