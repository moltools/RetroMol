"""Distinct assembly projections of an explored reaction hypergraph.

A frontier stops at an identified node or an unresolved leaf; otherwise it
chooses one reaction and includes ALL its products. Different derivations of
the same frontier count once. Atom assignments distinguish repeated monomers.
"""

import random
from functools import cached_property

import networkx as nx

from retromol.chem.tagging import get_tags_mol
from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.reaction_graph import ReactionGraph


Frontier = frozenset[str]


class AssemblyLimitError(RuntimeError):
    """Exact enumeration exceeded its budget; no partial count is returned."""


class AssemblyConstraintError(ValueError):
    """No explored assembly satisfies the terminal-unit constraints."""


class AssemblySpace:
    """Enumerate once, then count, inspect or sample distinct assemblies."""

    def __init__(
        self, reaction_graph: ReactionGraph, root_enc: str, *,
        max_assemblies: int = 10_000, max_combinations: int = 1_000_000,
        input_heavy_atoms: int | None = None,
    ):
        if root_enc not in reaction_graph.nodes:
            raise ValueError(f"Root {root_enc!r} is absent from the reaction graph")
        if max_assemblies < 1 or max_combinations < 1:
            raise ValueError("Enumeration limits must be positive")
        self.reaction_graph = reaction_graph
        self.root_enc = root_enc
        self.max_assemblies = max_assemblies
        self.max_combinations = max_combinations
        root_heavy_atoms = reaction_graph.nodes[root_enc].mol.GetNumHeavyAtoms()
        if input_heavy_atoms is None:
            input_heavy_atoms = root_heavy_atoms
        if type(input_heavy_atoms) is not int or input_heavy_atoms < root_heavy_atoms:
            raise ValueError("Input heavy-atom count must be an integer at least as large as the root's count")
        self.input_heavy_atoms = input_heavy_atoms

    @cached_property
    def _node_tags(self) -> dict[str, frozenset[int]]:
        root_tags = get_tags_mol(self.reaction_graph.nodes[self.root_enc].mol)
        return {enc: frozenset(get_tags_mol(node.mol) & root_tags)
                for enc, node in self.reaction_graph.nodes.items()}

    @cached_property
    def _terminal_tags(self) -> dict[str, frozenset[int]]:
        """Terminal regions reachable from this root, stopping at terminal nodes."""
        graph = self.reaction_graph
        pending = [self.root_enc]
        seen = set()
        terminal = {}
        while pending:
            enc = pending.pop()
            if enc in seen:
                continue
            seen.add(enc)
            node = graph.nodes[enc]
            if node.is_identified and node.identity is not None and node.identity.terminal:
                if self._node_tags[enc]:
                    terminal[enc] = self._node_tags[enc]
                continue
            for index in graph.out_edges.get(enc, []):
                edge = graph.edges[index]
                if enc not in edge.dsts:
                    pending.extend(edge.dsts)
        return terminal

    def preserves_terminals(self, frontier: Frontier) -> bool:
        """Whether a selection keeps terminal regions intact across all routes."""
        selected = [self._node_tags[enc] for enc in frontier]
        selected_terminals = [self._terminal_tags[enc] for enc in frontier if enc in self._terminal_tags]
        for region in set(self._terminal_tags.values()):
            if any(region <= tags for tags in selected):
                continue
            if any(region & tags and tags - region for tags in selected_terminals):
                continue
            return False
        return True

    @cached_property
    def frontiers(self) -> tuple[Frontier, ...]:
        """Distinct selections that preserve terminal units, at all coverages."""
        return tuple(frontier for frontier in self._raw_frontiers if self.preserves_terminals(frontier))

    @cached_property
    def _raw_frontiers(self) -> tuple[Frontier, ...]:
        """Enumerate routes before enforcing terminal boundaries across routes."""
        graph = self.reaction_graph
        tags = self._node_tags

        # Memoization only depends on ancestors in the same strongly connected
        # component. Ordinary acyclic branches share all their cached results.
        dependencies = nx.DiGraph()
        dependencies.add_nodes_from(graph.nodes)
        for edge in graph.edges:
            if edge.src not in edge.dsts:
                dependencies.add_edges_from((edge.src, child) for child in edge.dsts)
        components = {}
        for component in nx.strongly_connected_components(dependencies):
            members = frozenset(component)
            for enc in members:
                components[enc] = members

        memo: dict[tuple[str, frozenset[str]], dict[Frontier, frozenset[int]]] = {}
        combinations = 0

        def check_size(options: dict) -> None:
            if len(options) > self.max_assemblies:
                raise AssemblyLimitError(
                    f"Exact assembly enumeration exceeded max_assemblies={self.max_assemblies} "
                    "at an intermediate state; increase the limit to count or sample this graph."
                )

        def visit(enc: str, ancestors: frozenset[str]) -> dict[Frontier, frozenset[int]]:
            nonlocal combinations
            if enc in ancestors:
                return {}  # no cyclic derivations
            key = (enc, ancestors & components[enc])
            if key in memo:
                return memo[key]
            if not tags[enc]:
                return {frozenset(): frozenset()}  # reaction-introduced byproduct
            if len(ancestors) >= 400:
                raise AssemblyLimitError("Exact assembly enumeration exceeded the 400-step depth limit")
            node = graph.nodes[enc]
            stop = {frozenset([enc]): tags[enc]}
            options = dict(stop) if node.is_identified else {}
            terminal = node.is_identified and node.identity is not None and node.identity.terminal
            if not terminal:
                for index in graph.out_edges.get(enc, []):
                    edge = graph.edges[index]
                    if enc in edge.dsts or not edge.dsts:
                        continue
                    partial: dict[Frontier, frozenset[int]] = {frozenset(): frozenset()}
                    for child in edge.dsts:
                        child_options = visit(child, ancestors | {enc})
                        joined = {}
                        for frontier, occupied in partial.items():
                            for child_frontier, child_tags in child_options.items():
                                combinations += 1
                                if combinations > self.max_combinations:
                                    raise AssemblyLimitError(
                                        "Exact assembly enumeration exceeded "
                                        f"max_combinations={self.max_combinations}; increase the limit."
                                    )
                                if not occupied.isdisjoint(child_tags):
                                    continue
                                joined[frontier | child_frontier] = occupied | child_tags
                                check_size(joined)
                        partial = joined
                        if not partial:
                            break
                    options.update(partial)
                    check_size(options)
            # A fragment with no finite, compatible expansion remains unresolved.
            # This also preserves original atoms at dead ends of cyclic branches.
            if not options:
                options = stop
            memo[key] = options
            return options

        options = visit(self.root_enc, frozenset())
        return tuple(sorted(options, key=lambda frontier: tuple(sorted(frontier))))

    @cached_property
    def _coverages(self) -> tuple[float, ...]:
        """Score frontiers without constructing all their assembly graphs."""
        graph = self.reaction_graph
        root_heavy_tags = {
            atom.GetIsotope() for atom in graph.nodes[self.root_enc].mol.GetAtoms()
            if atom.GetAtomicNum() > 1 and atom.GetIsotope() != 0
        }
        identified_tags = {
            enc: get_tags_mol(node.mol) & root_heavy_tags
            for enc, node in graph.nodes.items() if node.is_identified
        }
        return tuple(
            len(set().union(*(identified_tags[enc] for enc in frontier if enc in identified_tags)))
            / self.input_heavy_atoms if self.input_heavy_atoms else 0.0
            for frontier in self.frontiers
        )

    @cached_property
    def max_coverage(self) -> float:
        """Highest attainable coverage among all compatible assembly views."""
        return max(self._coverages, default=0.0)

    def eligible_frontiers(self, min_coverage: float | None = None) -> tuple[Frontier, ...]:
        """Keep maximum-coverage views by default, or views meeting a minimum.

        An explicit minimum is an inclusive fraction from 0 to 1. Zero includes
        all views; an unattainable minimum returns no views. Filtering happens
        after exact enumeration, so the enumeration budgets still apply.
        """
        if min_coverage is not None and not 0 <= min_coverage <= 1:
            raise ValueError("Minimum coverage must be a finite fraction between 0 and 1")
        if min_coverage == 0:
            return self.frontiers
        threshold = self.max_coverage if min_coverage is None else min_coverage
        return tuple(frontier for frontier, coverage in zip(self.frontiers, self._coverages)
                     if coverage >= threshold)

    def count(self, *, min_coverage: float | None = 0.0) -> int:
        """Exact count: all views by default; None counts only maximum coverage."""
        return len(self.eligible_frontiers(min_coverage))

    def build(self, frontier: Frontier) -> AssemblyGraph:
        """Build a selected frontier using bonds of the stored root molecule."""
        if frontier not in self.frontiers:
            raise ValueError("Selection is not a valid assembly frontier")
        graph = self.reaction_graph
        return AssemblyGraph.build(
            root_mol=graph.nodes[self.root_enc].mol,
            monomers=[graph.nodes[enc] for enc in sorted(frontier)],
            include_unassigned=True,
        )

    def reaction_graph_for(self, frontier: Frontier) -> ReactionGraph:
        """Recover one explored route for a valid frontier, without inventing edges.

        Used when the cost-based default route violates terminal boundaries.
        Every branch must yield exactly its assigned frontier nodes; cyclic
        derivations and overlapping product assignments are rejected.
        """
        if frontier not in self.frontiers:
            raise ValueError("Selection is not a valid assembly frontier")
        graph = self.reaction_graph
        tags = self._node_tags
        required = {
            enc: frozenset(target for target in frontier if tags[target] & node_tags)
            for enc, node_tags in tags.items()
        }
        memo = {}
        attempts = 0

        def visit(enc: str, ancestors: frozenset[str]) -> frozenset[int] | None:
            nonlocal attempts
            if enc in ancestors:
                return None
            if len(ancestors) >= 400:
                raise AssemblyLimitError("Assembly route recovery exceeded the 400-step depth limit")
            key = (enc, ancestors)
            if key in memo:
                return memo[key]
            wanted = required[enc]
            if enc in frontier or not tags[enc]:
                return frozenset()
            if enc in self._terminal_tags or any(not tags[target] <= tags[enc] for target in wanted):
                return None
            for index in graph.out_edges.get(enc, []):
                attempts += 1
                if attempts > self.max_combinations:
                    raise AssemblyLimitError("Assembly route recovery exceeded max_combinations="
                                             f"{self.max_combinations}; increase the limit.")
                edge = graph.edges[index]
                if enc in edge.dsts or not edge.dsts:
                    continue
                assigned = frozenset()
                occupied = frozenset()
                edges = frozenset([index])
                for child in edge.dsts:
                    if occupied & tags[child] or assigned & required[child]:
                        break
                    child_edges = visit(child, ancestors | {enc})
                    if child_edges is None:
                        break
                    occupied |= tags[child]
                    assigned |= required[child]
                    edges |= child_edges
                else:
                    if assigned == wanted:
                        memo[key] = edges
                        return edges
            memo[key] = None
            return None

        indices = visit(self.root_enc, frozenset())
        if indices is None:
            raise AssemblyConstraintError("No explored reaction route yields the selected assembly")
        selected = ReactionGraph()
        nodes = {self.root_enc}
        for index in sorted(indices):
            edge = graph.edges[index]
            nodes.add(edge.src)
            nodes.update(edge.dsts)
            selected.out_edges.setdefault(edge.src, []).append(len(selected.edges))
            selected.edges.append(edge)
        for enc in sorted(nodes):
            selected.nodes[enc] = graph.nodes[enc]
            selected.out_edges.setdefault(enc, [])
        return selected

    def sample(
        self, n: int = 1, *, seed: int | None = None, min_coverage: float | None = None,
    ) -> list[AssemblyGraph]:
        """Uniform sampling among eligible views, without replacement."""
        if type(n) is not int or n < 0:
            raise ValueError("Sample size must be a nonnegative integer")
        if n == 0:
            return []
        eligible = self.eligible_frontiers(min_coverage)
        chosen = random.Random(seed).sample(eligible, min(n, len(eligible)))
        return [self.build(frontier) for frontier in chosen]


def count_unique_assemblies(
    reaction_graph: ReactionGraph, root_enc: str, *, min_coverage: float | None = 0.0, **limits,
) -> int:
    """Count all views, or qualifying views with the given coverage policy."""
    return AssemblySpace(reaction_graph, root_enc, **limits).count(min_coverage=min_coverage)


def sample_assemblies(
    reaction_graph: ReactionGraph, root_enc: str, n: int = 1, *, seed: int | None = None,
    min_coverage: float | None = None, **limits,
) -> list[AssemblyGraph]:
    """Sample distinct assemblies, defaulting to the highest coverage only."""
    return AssemblySpace(reaction_graph, root_enc, **limits).sample(n, seed=seed, min_coverage=min_coverage)
