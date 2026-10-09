"""Regression tests for extracting usable results from cyclic reaction graphs."""

import networkx as nx
import pytest

from retromol.chem.mol import encode_mol, smiles_to_mol
from retromol.model.reaction_graph import ReactionGraph, ReactionStep
from retromol.model.readout import LinearReadout
from retromol.model.result import Result
from retromol.model.rules import MatchingRule, RuleSet
from retromol.model.submission import Submission
from retromol.pipelines.parsing import extract_min_edge_synthesis_subgraph, run_retromol_with_timeout


def add_reaction(graph: ReactionGraph, reactant: str, *products: str) -> None:
    source = graph.add_node(smiles_to_mol(reactant))
    graph.add_edge(
        source,
        [smiles_to_mol(smi) for smi in products],
        ReactionStep(kind="contested", names=("test reaction",), rule_ids=("test",)),
    )


def identify(graph: ReactionGraph, smiles: str) -> None:
    node = graph.nodes[encode_mol(smiles_to_mol(smiles))]
    node.identify([MatchingRule(name="known", smiles=smiles, props={}, pseudonyms=[])])


def assert_acyclic(graph: ReactionGraph) -> None:
    directed = nx.DiGraph()
    directed.add_nodes_from(graph.nodes)
    directed.add_edges_from((edge.src, dst) for edge in graph.edges for dst in edge.dsts)
    assert nx.is_directed_acyclic_graph(directed)


def test_self_dependent_reaction_retains_unresolved_root() -> None:
    graph = ReactionGraph()
    add_reaction(graph, "[1CH4]", "[1CH4]", "O")
    identify(graph, "O")
    root = encode_mol(smiles_to_mol("[1CH4]"))

    result = extract_min_edge_synthesis_subgraph(graph, root)

    assert set(result.graph.nodes) == {root}
    assert result.graph.edges == []
    assert not result.solved
    assert result.total_cost == 5.0


def test_cyclic_fragment_does_not_discard_identified_siblings() -> None:
    submission = Submission("CCO")
    graph = ReactionGraph()
    add_reaction(graph, "[1CH3][2CH2][3OH]", "[1CH4]", "[2CH3][3OH]")
    add_reaction(graph, "[2CH3][3OH]", "[2CH2]=[3O]")
    add_reaction(graph, "[2CH2]=[3O]", "[2CH3][3OH]")
    identify(graph, "[1CH4]")
    root = encode_mol(submission.mol)

    extracted = extract_min_edge_synthesis_subgraph(graph, root)
    result = Result(submission, extracted.graph, LinearReadout.from_reaction_graph(root, extracted.graph))

    assert root in extracted.graph.nodes
    assert not extracted.solved
    assert_acyclic(extracted.graph)
    assert len(extracted.graph.get_leaf_nodes(identified_only=True)) == 1
    assert result.calculate_coverage() == pytest.approx(1 / 3)


def test_cycle_with_identified_exit_is_solved() -> None:
    graph = ReactionGraph()
    add_reaction(graph, "[1CH4]", "[1CH3]O")
    add_reaction(graph, "[1CH3]O", "[1CH4]")
    add_reaction(graph, "[1CH3]O", "[1CH3]N")
    identify(graph, "[1CH3]N")
    root = encode_mol(smiles_to_mol("[1CH4]"))

    result = extract_min_edge_synthesis_subgraph(graph, root)

    assert result.solved
    assert len(result.graph.edges) == 2
    assert_acyclic(result.graph)


def test_fully_identified_route_can_cost_more_than_unsolved_penalty() -> None:
    graph = ReactionGraph()
    add_reaction(graph, "[1CH4]", "[1CH3]O")
    add_reaction(graph, "[1CH3]O", "[1CH3]N")
    add_reaction(graph, "[1CH3]N", "[1CH3]Cl")
    identify(graph, "[1CH3]Cl")
    root = encode_mol(smiles_to_mol("[1CH4]"))

    result = extract_min_edge_synthesis_subgraph(graph, root, edge_base_cost=2.0)

    assert result.total_cost > 5.0
    assert result.solved


@pytest.mark.parametrize(
    "smiles",
    [
        "CCOS(=O)(=O)O",
        # NPA000139 and NPA000469 reproduced the missing-root error with default rules.
        "O=C(CC(O)(CC(=O)NCCCNCCCNC(=O)c1ccc(O)c(O)c1)C(=O)O)NCCCCNCCCNC(=O)c1ccc(O)c(O)c1S(=O)(=O)O",
        "CCCCCCCCCCCCCCCC(=O)NC1CC(O)C(O)NC(=O)C2C(O)C(C)CN2C(=O)C(C(O)CC(N)=O)NC(=O)C(C(O)Cc2ccc(OS(=O)(=O)O)c(O)c2)NC(=O)C2CC(O)CN2C(=O)C(C(C)O)NC1=O",
    ],
    ids=["ethyl sulfate", "NPA000139", "NPA000469"],
)
def test_default_rules_preserve_partial_results_for_sulfated_compounds(smiles: str, ruleset: RuleSet) -> None:
    result = run_retromol_with_timeout(Submission(smiles), ruleset)

    assert encode_mol(result.submission.mol) in result.reaction_graph.nodes
    assert_acyclic(result.reaction_graph)
    assert 0 < result.calculate_coverage() < 1
    assert any(not node.is_identified for node in result.linear_readout.assembly_graph.monomer_nodes())
    restored = Result.from_dict(result.to_dict())
    assert restored.calculate_coverage() == result.calculate_coverage()


def test_missing_root_in_invalid_input_graph_still_raises() -> None:
    with pytest.raises(ValueError, match="Root encoding missing"):
        extract_min_edge_synthesis_subgraph(ReactionGraph(), "missing")
