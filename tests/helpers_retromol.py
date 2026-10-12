"""Shared helpers for RetroMol integration tests."""

from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import networkx as nx
import yaml

from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.rules import RuleSet
from retromol.model.result import Result
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol_with_timeout


def load_rule_set() -> RuleSet:
    """
    Load the default RetroMol rule set once.

    :return: the loaded RuleSet object
    """
    return RuleSet.load_default(match_stereochemistry=False)


def load_stereo_rule_set() -> RuleSet:
    """
    Load the default RetroMol rule set once, matching stereochemistry.

    :return: the loaded RuleSet object
    """
    return RuleSet.load_default(match_stereochemistry=True)


def parse_compound(smiles: str, ruleset: RuleSet) -> Result:
    """
    Parse a compound SMILES string into an Result object.

    :param smiles: the SMILES string of the compound to parse
    :param ruleset: the RuleSet to use for parsing
    :return: the resulting Result object
    """
    submission = Submission(smiles)
    return run_retromol_with_timeout(submission, ruleset)


@dataclass(frozen=True)
class ExpectedAssembly:
    """An explicit, node-ID-independent assembly graph fixture.

    Nodes are (monomer name or None, original heavy-atom count). None denotes
    an unresolved region. Edges are (node index, node index, root bond labels),
    e.g. (0, 1, ("C-N:SINGLE",)). Repeated bond labels represent multiple root
    bonds between the same two monomers, as in a cyclic dipeptide.
    """

    coverage: float
    nodes: tuple[tuple[str | None, int], ...]
    edges: tuple[tuple[int, int, tuple[str, ...]], ...]

    def graph(self) -> nx.Graph:
        if not 0 <= self.coverage <= 1:
            raise ValueError("Expected coverage must be between 0 and 1")
        graph = nx.Graph()
        for index, (name, heavy_atoms) in enumerate(self.nodes):
            if (name is not None and not isinstance(name, str)) or heavy_atoms < 1:
                raise ValueError(f"Invalid expected assembly node: {(name, heavy_atoms)!r}")
            graph.add_node(index, name=name, heavy_atoms=heavy_atoms)
        for left, right, bonds in self.edges:
            if left not in graph or right not in graph or left == right or not bonds:
                raise ValueError(f"Invalid expected assembly edge: {(left, right, bonds)!r}")
            if graph.has_edge(left, right):
                raise ValueError("Put all bonds between two monomers in a single edge")
            graph.add_edge(left, right, bonds=tuple(sorted(bonds)))
        return graph


def load_expected_assemblies(path: Path) -> dict[str, tuple[ExpectedAssembly, ...]]:
    """Read fixed, reviewable fixtures; never derive expectations during a test."""
    raw = yaml.safe_load(path.read_text())
    expectations = {}
    for name, alternatives in raw.items():
        if not alternatives:
            raise ValueError(f"No expected assemblies supplied for {name!r}")
        expected = tuple(ExpectedAssembly(
            coverage=alternative["coverage"],
            nodes=tuple((label, count) for label, count in alternative["nodes"]),
            edges=tuple((left, right, tuple(bonds)) for left, right, bonds in alternative["edges"]),
        ) for alternative in alternatives)
        for alternative in expected:
            alternative.graph()  # Reject malformed fixtures during collection.
        expectations[name] = expected
    return expectations


def assembly_comparison_graph(result: Result, assembly: AssemblyGraph) -> nx.Graph:
    """Keep labels, original atom counts and connectivity, excluding empty nodes.

    Node hashes, root atom numbering, reaction history, and reaction-added atoms
    are not assembly identities. Nonempty unresolved regions and every root
    bond between regions remain part of the comparison.
    """
    root_heavy_tags = {
        atom.GetIsotope() for atom in result.root.mol.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetIsotope() != 0
    }
    graph = nx.Graph()
    for key, data in assembly.g.nodes(data=True):
        heavy_atoms = len(data["tags"] & root_heavy_tags)
        if not heavy_atoms:
            continue
        node = data["molnode"]
        name = node.identity.name if node is not None and node.is_identified else None
        graph.add_node(key, name=name, heavy_atoms=heavy_atoms)
    for left, right, links in assembly.edges_with_bonds():
        if left in graph and right in graph:
            bonds = tuple(sorted(
                f'{"-".join(sorted((link.a1_symbol, link.a2_symbol)))}:{link.bond_type}'
                for link in links
            ))
            graph.add_edge(left, right, bonds=bonds)
    return graph


def _graph_summary(graph: nx.Graph) -> str:
    """Include numbered edges so equal monomer lists with wrong wiring are clear."""
    indices = {node: index for index, node in enumerate(graph)}
    nodes = [(data["name"], data["heavy_atoms"]) for _, data in graph.nodes(data=True)]
    edges = [(indices[left], indices[right], data["bonds"])
             for left, right, data in graph.edges(data=True)]
    return f"nodes={nodes!r}, edges={edges!r}"


def assert_assemblies(
    result: Result,
    expected: Sequence[ExpectedAssembly],
    *,
    coverage_tolerance: float = 1e-2,
) -> None:
    """Require ANY expected graph to occur among ALL valid parsed assemblies.

    Extra assemblies are allowed. Search is exhaustive and deterministic, not
    random sampling or selection of the default route. Enumeration-limit errors
    propagate instead of silently accepting an incomplete search.
    """
    if not expected:
        raise ValueError("At least one expected assembly graph is required")
    expected_graphs = [(alternative, alternative.graph()) for alternative in expected]
    space = result.assembly_space()
    frontiers = space.eligible_frontiers(
        min_coverage=max(0.0, min(alternative.coverage for alternative in expected) - coverage_tolerance),
    )
    diagnostics = []
    node_match = nx.algorithms.isomorphism.categorical_node_match(
        ["name", "heavy_atoms"], [None, None],
    )
    edge_match = nx.algorithms.isomorphism.categorical_edge_match("bonds", ())
    for frontier in frontiers:
        assembly = space.build(frontier)
        coverage = result.calculate_coverage(assembly)
        actual_graph = assembly_comparison_graph(result, assembly)
        for alternative, expected_graph in expected_graphs:
            if abs(coverage - alternative.coverage) > coverage_tolerance:
                continue
            if nx.is_isomorphic(actual_graph, expected_graph, node_match=node_match, edge_match=edge_match):
                return
        if len(diagnostics) < 5:
            diagnostics.append(f"coverage={coverage:.6f}, {_graph_summary(actual_graph)}")

    expected_descriptions = [
        f"coverage={alternative.coverage:.6f}, {_graph_summary(graph)}"
        for alternative, graph in expected_graphs
    ]
    raise AssertionError(
        "No parsed assembly matches any expected graph.\nExpected alternatives:\n  "
        + "\n  ".join(expected_descriptions)
        + f"\nParsed {len(space.frontiers)} valid assemblies; {len(frontiers)} meet the coverage lower bound."
        + f" Maximum coverage: {space.max_coverage:.6f}.\nFirst candidates (at most 5):\n  "
        + "\n  ".join(diagnostics)
    )
