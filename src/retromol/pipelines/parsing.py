"""Module for applying reaction rules to molecules using a reaction graph approach."""

import logging
import os
from math import inf
from collections import deque, defaultdict
from typing import Dict, List, Set, Optional

from retromol.utils.timeout import timeout_decorator
from retromol.model.submission import Submission
from retromol.model.rules import RuleSet, index_uncontested, apply_uncontested
from retromol.model.result import Result
from retromol.model.reaction_graph import ReactionGraph, ReactionStep, RxnEdge
from retromol.model.synthesis import SynthesisExtractResult
from retromol.chem.mol import Mol, encode_mol
from retromol.chem.tagging import get_tags_mol
from retromol.visualization.reaction_graph import visualize_reaction_graph


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

        # If we've already expanded this encoding, skip
        if parent_enc in expanded:
            continue
        expanded.add(parent_enc)

        # Identity gating
        node = g.nodes[parent_enc]
        ident = node.identify(matching_rules, match_stereochemistry=False)
        if ident:
            continue

        # Uncontested in bulk (combined step)
        uncontested = index_uncontested(parent, reaction_rules, failed_combos)
        if uncontested:
            products, applied_in_bulk, new_failed = apply_uncontested(parent, uncontested, original_taken_tags)
            failed_combos.update(new_failed)

            # If uncontested existed but none succeed, fall through to contested
            if applied_in_bulk:
                step = ReactionStep(
                    kind="uncontested",
                    names=tuple(rl.name for rl in applied_in_bulk),
                    rule_ids=tuple(rl.id for rl in applied_in_bulk),
                )
                g.add_edge(parent_enc, products, step)

                # Enqueue newly discovered products (by encoding)
                for m in products:
                    enc = encode_mol(m)
                    if enc not in expanded and enc not in enqueued:
                        q.append(Mol(m))
                        enqueued.add(enc)

                continue  # skip contested for this parent if bulk succeeded

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
    root_enc: int,
    prefer_kind: tuple[str, ...] = ("uncontested", "contested"),
    edge_base_cost: float = 1.0,
) -> SynthesisExtractResult:
    """
    Extract a minimum-edge synthesis subgraph from a retrosynthesis ReactionGraph.

    Interprets the graph as an AND/OR graph:
      - molecule node: OR (choose one outgoing reaction edge)
      - reaction edge: AND (must solve all dst precursor nodes)
      - identified molecule nodes are terminal solved leaves (cost=0)

    The extracted subgraph contains at most one chosen outgoing edge per expanded node
    and includes all required precursor branches for that choice.
    """
    # Adjacency list of outgoing edges for quick access
    out_edges: Dict[int, List[int]] = defaultdict(list)
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
    memo_cost: Dict[int, float] = {}
    memo_choice: Dict[int, Optional[int]] = {}  # node -> chosen edge index
    visiting: Set[int] = set()

    def is_identified(enc: int) -> bool:
        n = g.nodes.get(enc)
        return bool(n and n.is_identified)

    def solve_cost(enc: int) -> float:
        # Identified nodes are leaves
        if is_identified(enc):
            memo_choice[enc] = None
            return 0.0

        if enc in memo_cost:
            return memo_cost[enc]

        if enc in visiting:
            # cycle guard: treat as unsolvable in this simple DP
            return inf

        visiting.add(enc)

        best = inf
        best_ei: Optional[int] = None

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
    if total == inf:
        return SynthesisExtractResult(graph=ReactionGraph(), solved=False, total_cost=inf)

    # Extract chosen policy edges into a new small graph
    new_g = ReactionGraph()

    kept_nodes: Set[int] = set()
    kept_edge_indices: Set[int] = set()

    def extract(enc: int) -> None:
        if enc in kept_nodes:
            return
        kept_nodes.add(enc)

        # Always keep the node
        if enc in g.nodes:
            new_g.nodes[enc] = g.nodes[enc]
            new_g.out_edges.setdefault(enc, [])

        # Stop at identified leaves
        if is_identified(enc):
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

    return SynthesisExtractResult(graph=new_g, solved=True, total_cost=total)


def run_retromol(submission: Submission, rules: RuleSet) -> Result:
    """
    Run RetroMol retrosynthesis on the given input molecule using the specified reaction rules.
    
    :param submission: Submission object containing the input molecule and data
    :param rules: Rules object containing the reaction rules to apply
    :return: Result object containing the retrosynthesis results
    """

    g = process_mol(submission, rules)
    print(g)

    root = encode_mol(submission.mol)
    r = extract_min_edge_synthesis_subgraph(g, root_enc=root)
    print(r.graph)

    # store in downloads
    visualize_reaction_graph(r.graph, html_path="/Users/davidmeijer/Downloads/reaction_graph.html")

    return Result()


run_retromol_with_timeout = timeout_decorator(seconds=int(os.getenv("TIMEOUT_RUN_RETROMOL", "30")))(run_retromol)
