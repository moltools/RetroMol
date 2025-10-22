"""This module provides the main API for the RetroMol package."""

import logging
from importlib.resources import files
from typing import Any, cast

import yaml
from networkx import Graph

import retromol.data
from retromol.apply import resolve_mol
from retromol.chem import get_tags_mol, smiles_to_mol
from retromol.config import LOGGER_NAME, TIMEOUT_RUN_RETROMOL
from retromol.helpers import timeout_decorator
from retromol.io import Input, Result
from retromol.rules import Rules, load_rules_from_files


def find_eligible_nodes(
    graph: "Graph[int | str]", parse_identified_nodes: bool = False
) -> list[tuple["Graph[int | str]", int | str]]:
    """
    Walk `graph` and nested sub-graphs (via node attr 'graph') and return a list of
    (parent_graph, node_id) for nodes whose attrs['graph'] is None, but ONLY when
    their parent graph is at the current deepest nesting level.

    Effectively: expand the most recently created level of graphs, one wave at a time.

    If parse_identified_nodes is False (default), nodes with a non-None 'identity'
    are skipped; if True, identity is ignored.

    :param graph: the root graph to search
    :param parse_identified_nodes: whether to include nodes with an 'identity' attribute
    :return: a list of tuples containing the parent graph and node ID of eligible nodes
    """

    # Determine the maximum depth of any existing graph in the nesting.
    # Depth of the root `graph` is 0; a subgraph inside a node of depth d has depth d+1.
    def _max_depth(g: "Graph[int | str]", depth: int = 0) -> int:
        md = depth
        for _, attrs in g.nodes(data=True):
            sub: Any = attrs.get("graph")
            if isinstance(sub, Graph):
                sub = cast("Graph[int | str]", sub)
                md = max(md, _max_depth(sub, depth + 1))
        return md

    max_depth = _max_depth(graph)

    # Collect eligible nodes whose parent graph is exactly at max_depth.
    eligible: list[tuple[Graph[int | str], int | str]] = []

    def _collect(g: "Graph[int | str]", depth: int = 0):
        # Only collect from graphs at the current frontier depth.
        if depth == max_depth:
            for nid, attrs in g.nodes(data=True):
                if attrs.get("graph") is None:
                    if parse_identified_nodes or attrs.get("identity") is None:
                        eligible.append((g, nid))
        # Recurse to find deeper graphs (to maintain correctness of traversal),
        # but we only collect at exactly max_depth.
        for _, attrs in g.nodes(data=True):
            sub = attrs.get("graph")
            if isinstance(sub, Graph):
                sub = cast("Graph[int | str]", sub)
                _collect(sub, depth + 1)

    _collect(graph)
    return eligible


def run_retromol(
    input: Input,
    rule_set: Rules | None = None,
    wave_configs: list[dict[str, Any]] | None = None,
    match_stereochemistry: bool = False,
) -> Result:
    """
    Run RetroMol on a given input compound.

    :param input: the input compound to process, as an Input object
    :param rule_set: the set of rules to apply, as a Rules object
    :param wave_configs: configuration for each wave, as a dictionary mapping wave numbers to config dicts
    :param match_stereochemistry: whether to match stereochemistry during processing
    :return: a Result object containing the processed graph and metadata
    """
    if rule_set is None:
        path_rx = str(files(retromol.data).joinpath("default_reaction_rules.yml"))
        path_mx = str(files(retromol.data).joinpath("default_matching_rules.yml"))
        rule_set = load_rules_from_files(path_rx, path_mx)

    if wave_configs is None:
        path_wave_config = str(files(retromol.data).joinpath("default_wave_config.yml"))
        with open(path_wave_config) as f:
            wave_configs = yaml.safe_load(f)

    logger = logging.getLogger(LOGGER_NAME)
    matching_rules = rule_set.get_matching_rules()

    logger.debug(f"Processing input: {input.cid} with SMILES: {input.smi}")

    props = input.props if input.props else {}
    motif_graph = None

    reserved_tags: set[int] = set(get_tags_mol(input.mol))

    # Loop through wave numbers and apply rules to the motif graph, increasing the motif nesting
    for wave_config in wave_configs or []:
        if motif_graph is None:
            # First wave: create the initial motif graph from the input molecule
            logger.debug(f"Starting wave {wave_config.get('wave_name', 'unnamed')} on input molecule")
            motif_graph = resolve_mol(
                input.mol,
                reserved_tags,
                rule_set.get_reaction_rules(group_names=wave_config.get("reaction_groups", [])),
                matching_rules,
                match_stereochemistry=match_stereochemistry,
                wave_config=wave_config,
            )
            continue

        # Find all nodes that are eligible for processing: all nodes without an identity and sub-graph
        parse_identified_nodes = wave_config.get("parse_identified_nodes", False)
        todo = find_eligible_nodes(motif_graph, parse_identified_nodes=parse_identified_nodes)

        # Loop through the eligible nodes and process them
        for parent_graph, n_id in todo:
            attrs = parent_graph.nodes[n_id]
            logger.debug(f"Processing node {n_id} in graph with {parent_graph.number_of_nodes()} nodes")
            sub = resolve_mol(
                smiles_to_mol(attrs["smiles"]),
                reserved_tags,
                rule_set.get_reaction_rules(group_names=wave_config.get("reaction_groups", [])),
                matching_rules,
                match_stereochemistry=match_stereochemistry,
                wave_config=wave_config,
            )
            attrs["graph"] = sub

    # Return the motif graph as a labeled result
    return Result(
        input_id=input.cid,
        graph=motif_graph if motif_graph is not None else Graph(),
        props=props,
        sha256_reaction_rules=rule_set.sha256_reaction_rules,
        sha256_matching_rules=rule_set.sha256_matching_rules,
    )


run_retromol_with_timeout = timeout_decorator(seconds=TIMEOUT_RUN_RETROMOL)(run_retromol)
