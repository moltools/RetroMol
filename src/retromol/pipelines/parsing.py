"""Module for applying reaction rules to molecules using a reaction graph approach."""

import logging
import os
from math import inf
from collections import deque, defaultdict
from typing import Optional

from retromol.utils.timeout import timeout_decorator
from retromol.model.submission import Submission
from retromol.model.rules import RuleSet, index_uncontested, apply_uncontested
from retromol.model.result import Result
from retromol.model.reaction_graph import ReactionGraph, ReactionStep, RxnEdge
from retromol.model.readout import LinearReadout
from retromol.model.synthesis import SynthesisExtractResult
from retromol.chem.mol import Mol, encode_mol, mol_to_smiles
from retromol.chem.tagging import get_tags_mol


log = logging.getLogger(__name__)


def process_mol(submission: Submission, ruleset: RuleSet) -> ReactionGraph:
    """
    Process a molecule by applying reaction rules and constructing a reaction graph.

    :param submission: Submission: the input molecule and associated data
    :param ruleset: RuleSet: the set of reaction and matching rules to apply
    :return: ReactionGraph: the resulting reaction graph after processing
    """
    reaction_rules = ruleset.reaction_rules
    matching_rules = ruleset.matching_rules

    g = ReactionGraph()

    original_taken_tags = get_tags_mol(submission.mol)
    failed_combos: set[tuple[int, frozenset[int]]] = set()

    # Track queue/expansion status by encoding to avoid duplicate work
    enqueued: set[str] = set()
    expanded: set[str] = set()

    q: deque[Mol] = deque()
    q.append(Mol(submission.mol))
    enqueued.add(encode_mol(submission.mol))

    while q:
        parent = q.popleft()
        parent_enc = g.add_node(parent)

        log.debug(f"expanding node {mol_to_smiles(parent, include_tags=False)}")

        # If we've already expanded this encoding, skip
        if parent_enc in expanded:
            continue
        expanded.add(parent_enc)

        # Identity gating: only gate (stop expanding) if identified AND terminal=True
        node = g.nodes[parent_enc]
        ident = node.identify(matching_rules, match_stereochemistry=ruleset.match_stereochemistry)
        if ident and bool(getattr(ident, "terminal", True)):
            log.debug(f"node (identity={ident.name}) identified as terminal; stopping expansion")
            continue

        # Uncontested in bulk (combined step)
        allowed_in_bulk = [rl for rl in reaction_rules if rl.allowed_in_bulk]
        uncontested = index_uncontested(parent, allowed_in_bulk, failed_combos)
        if uncontested:
            log.debug(f"applying {len(uncontested)} uncontested rule(s) in bulk")

            products, applied_in_bulk, new_failed = apply_uncontested(parent, uncontested, original_taken_tags)
            failed_combos.update(new_failed)

            # If uncontested existed but none succeed, fall through to contested
            if applied_in_bulk:
                step = ReactionStep(
                    kind="uncontested",
                    names=tuple(rl.name for rl, _ in applied_in_bulk),
                    rule_ids=tuple(rl.id for rl, _ in applied_in_bulk),
                )
                g.add_edge(parent_enc, products, step)

                # Enqueue newly discovered products (by encoding)
                for m in products:
                    enc = encode_mol(m)
                    if enc not in expanded and enc not in enqueued:
                        q.append(Mol(m))
                        enqueued.add(enc)

                continue

        # Contested exhaustive
        for rl in reaction_rules:

            results = rl.apply(parent, None)
            if not results:
                continue

            for result_set in results:
                # result_set is an iterable of product mols
                step = ReactionStep(
                    kind="contested",
                    names=(rl.name,),
                    rule_ids=(rl.id,),
                )
                g.add_edge(parent_enc, result_set, step)

                for m in result_set:
                    enc = encode_mol(m)
                    if enc not in expanded and enc not in enqueued:
                        q.append(Mol(m))
                        enqueued.add(enc)

    return g


def extract_min_edge_synthesis_subgraph(
    g: ReactionGraph,
    root_enc: str,
    prefer_kind: tuple[str, ...] = ("uncontested", "contested"),
    edge_base_cost: float = 1.0,
    nonterminal_leaf_penalty: float = 0.25,
    unsolved_leaf_penalty: float = 5.0,
) -> SynthesisExtractResult:
    """
    Extract a minimum-edge synthesis subgraph from a retrosynthesis ReactionGraph.

    Interprets the graph as an AND/OR graph:
      - molecule node: OR (choose one outgoing reaction edge)
      - reaction edge: AND (must solve all dst precursor nodes)
      - identified molecule nodes are terminal solved leaves (cost=0)

    The extracted subgraph contains at most one chosen outgoing edge per expanded node
    and includes all required precursor branches for that choice.

    :param g: ReactionGraph: the full retrosynthesis reaction graph
    :param root_enc: encoding of the root molecule to extract from
    :param prefer_kind: tuple[str, ...]: preference order for reaction kinds when costs are equal
    :param edge_base_cost: float: base cost per reaction edge
    :param nonterminal_leaf_penalty: float: penalty for identified leaves that are non-terminal
    :param unsolved_leaf_penalty: float: penalty for unsolved leaves (i.e., "give up" cost)
    :return: SynthesisExtractResult: the extracted synthesis subgraph and status
    :raises: ValueError: if root_enc is not in the graph
    """
    if root_enc not in g.nodes:
        raise ValueError(f"root encoding {root_enc} not found in reaction graph nodes")

    # Adjacency list of outgoing edges for quick access
    out_edges: dict[int, list[int]] = defaultdict(list)
    for ei, e in enumerate(g.edges):
        out_edges[e.src].append(ei)

    kind_rank = {k: i for i, k in enumerate(prefer_kind)}

    def edge_cost(e: RxnEdge) -> float:
        # Primary objective: fewer edges
        # Secondary tiebreakers: prefer uncontested; optionally penalize many precursors slightly
        kind_penalty = 0.001 * kind_rank.get(e.step.kind, 999)
        branch_penalty = 0.0001 * len(e.dsts)
        return edge_base_cost + kind_penalty + branch_penalty

    # DP memo: cost to "solve" a node into identified leaves
    memo_cost: dict[int, float] = {}
    memo_choice: dict[int, Optional[int]] = {}  # node -> chosen edge index
    visiting: set[int] = set()

    def is_terminal(enc: str) -> bool:
        n = g.nodes.get(enc)
        if not n or not n.is_identified:
            return False
        ident = n.identity
        return bool(getattr(ident, "terminal", True))

    def solve_cost(enc: str) -> float:
        n = g.nodes.get(enc)

        # Hard leaves: identified + terminal=True
        if n and n.is_identified and is_terminal(enc):
            memo_choice[enc] = None
            return 0.0
        
        # If nothing to expand, treat as "frontier leaf" (unsolved remainder)
        if not out_edges.get(enc):
            memo_choice[enc] = None
            # Identified nonterminal with no edges is still fine as 0, else unsolved penalty
            if n and n.is_identified:
                return 0.0
            return unsolved_leaf_penalty

        # Soft leaves: identified + terminal=False
        leaf_cost = inf
        if n and n.is_identified and not is_terminal(enc):
            # If no outgoing edges exist, must stop here (fallback leaf)
            if not out_edges.get(enc):
                memo_choice[enc] = None
                return 0.0
            # Otherwise, allow stopping, but discourage it
            leaf_cost = nonterminal_leaf_penalty

        if enc in memo_cost:
            return memo_cost[enc]

        if enc in visiting:
            # Cycle guard: treat as unsolvable in this simple DP
            return inf

        visiting.add(enc)

        best = leaf_cost
        best_ei: Optional[int] = None if best < inf else None

        for ei in out_edges.get(enc, []):
            e = g.edges[ei]
            # AND: all children must be solvable
            c = edge_cost(e)
            for d in e.dsts:
                dc = solve_cost(d)
                if dc == inf:
                    c = inf
                    break
                c += dc

            if c < best:
                best = c
                best_ei = ei

        visiting.remove(enc)

        memo_cost[enc] = best
        memo_choice[enc] = best_ei
        return best

    total = solve_cost(root_enc)
    # solved = no unsolved frontier was needed
    solved = total < inf and total < unsolved_leaf_penalty
    if total == inf:
        return SynthesisExtractResult(graph=ReactionGraph(), solved=False, total_cost=inf)

    # Extract chosen policy edges into a new small graph
    new_g = ReactionGraph()

    kept_nodes: set[int] = set()
    kept_edge_indices: set[int] = set()

    def extract(enc: str) -> None:
        if enc in kept_nodes:
            return
        kept_nodes.add(enc)

        # Always keep the node
        if enc in g.nodes:
            new_g.nodes[enc] = g.nodes[enc]
            new_g.out_edges.setdefault(enc, [])

        # Stop at identified leaves
        if is_terminal(enc):
            return

        ei = memo_choice.get(enc)
        if ei is None:
            return

        kept_edge_indices.add(ei)
        e = g.edges[ei]

        # Keep all dsts (AND)
        for d in e.dsts:
            extract(d)

    extract(root_enc)

    # Add edges (after nodes exist)
    for ei in kept_edge_indices:
        e = g.edges[ei]
        if e.src not in new_g.nodes:
            continue
        dsts = tuple(d for d in e.dsts if d in new_g.nodes)
        if not dsts:
            continue
        new_edge = RxnEdge(src=e.src, dsts=dsts, step=e.step)
        new_g.edges.append(new_edge)
        new_g.out_edges.setdefault(e.src, []).append(len(new_g.edges) - 1)

    return SynthesisExtractResult(graph=new_g, solved=solved, total_cost=total)


def run_retromol(submission: Submission, rules: RuleSet) -> Result:
    """
    Run RetroMol retrosynthesis on the given input molecule using the specified reaction rules.
    
    :param submission: Submission object containing the input molecule and data
    :param rules: Rules object containing the reaction rules to apply
    :return: Result object containing the retrosynthesis results
    """ 
    # Parse compound into reaction graph
    g = process_mol(submission, rules)
    log.debug(f"retrosynthesis graph has {len(g.nodes)} ({len(g.identified_nodes)} identified) nodes and {len(g.edges)} edges")

    # Extract minimum-edge synthesis subgraph
    root = encode_mol(submission.mol)
    r = extract_min_edge_synthesis_subgraph(
        g,
        root_enc=root,
        edge_base_cost=0.25,                # low base cost encourages longer syntheses
        nonterminal_leaf_penalty=100.0,     # high penalty forces expansion of non-terminal leaves (e.g., fatty acids)
    )
    log.debug(f"extracted synthesis subgraph has {len(r.graph.nodes)} ({len(r.graph.identified_nodes)} identified) nodes and {len(r.graph.edges)} edges")

    if not r.solved:
        log.debug("retrosynthesis extraction failed to find a solution")

    # Calculate the linear readouts for the synthesis graph
    linear_readout = LinearReadout.from_reaction_graph(root, r.graph)

    return Result(
        submission=submission,
        reaction_graph=r.graph,
        linear_readout=linear_readout,
    )


run_retromol_with_timeout = timeout_decorator(seconds=int(os.getenv("TIMEOUT_RUN_RETROMOL", "60")))(run_retromol)
