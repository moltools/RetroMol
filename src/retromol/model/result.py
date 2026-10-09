"""Compact parsing results, with assembly projections computed on demand."""

from dataclasses import dataclass
from functools import cached_property
from typing import Any

from retromol.chem.mol import Mol, encode_mol, mol_to_inchikey, smiles_to_mol
from retromol.chem.stereo import capture_double_bond_stereo
from retromol.chem.tagging import get_tags_mol, remove_tags
from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.assembly_sampling import AssemblyConstraintError, AssemblySpace
from retromol.model.submission import Submission
from retromol.model.reaction_graph import MolNode, ReactionGraph
from retromol.model.readout import LinearReadout


@dataclass(frozen=True)
class Result:
    """An explored reaction graph and its input/provenance metadata.

    Only these source data are serialized. Assemblies, coverage and sequences
    are projections requested by the consumer. The default projection uses a
    minimum-cost, fine-grained policy, weighting collapsed unresolved regions
    by their original heavy atoms.
    """

    submission: Submission
    reaction_graph: ReactionGraph
    reaction_rules_hash: str | None = None
    matching_rules_hash: str | None = None
    match_stereochemistry: bool = False
    root_enc: str | None = None

    def __post_init__(self) -> None:
        root = self.root_enc if self.root_enc is not None else encode_mol(self.submission.mol)
        if root not in self.reaction_graph.nodes:
            raise ValueError(f"Root {root!r} is absent from the reaction graph")
        object.__setattr__(self, "root_enc", root)

    def __str__(self) -> str:
        return f"Result(submission={self.submission}, reaction_graph={self.reaction_graph})"

    @property
    def root(self) -> MolNode:
        return self.reaction_graph.nodes[self.root_enc]

    @property
    def props(self) -> dict[str, Any] | None:
        return self.submission.props

    @cached_property
    def selected_reaction_graph(self) -> ReactionGraph:
        """One fine-grained route; selection never alters the stored graph."""
        from retromol.pipelines.parsing import extract_min_edge_synthesis_subgraph

        selected = extract_min_edge_synthesis_subgraph(
            self.reaction_graph, self.root_enc,
            edge_base_cost=0.25, nonterminal_leaf_penalty=100.0,
            weight_pruned_atoms=True,
        ).graph
        space = self.assembly_space()
        frontier = frozenset(node.enc for node in selected.get_leaf_nodes(identified_only=False))
        if space.preserves_terminals(frontier):
            return selected
        # A cheap route can cut a terminal unit before recognizing it. Choose
        # a valid, maximally covered alternative, preferring finer resolution,
        # and retain its actual reaction route for downstream reconstruction.
        candidates = space.eligible_frontiers()
        if not candidates:
            raise AssemblyConstraintError("No explored assembly preserves the identified terminal units")
        return space.reaction_graph_for(max(candidates, key=len))

    @cached_property
    def assembly_graph(self) -> AssemblyGraph:
        """Default assembly, calculated only when requested."""
        return AssemblyGraph.build(
            root_mol=self.root.mol,
            monomers=self.selected_reaction_graph.get_leaf_nodes(identified_only=False),
            include_unassigned=True,
        )

    @cached_property
    def linear_readout(self) -> LinearReadout:
        """Default sequence view; neither it nor its assembly is saved."""
        return LinearReadout.from_assembly_graph(self.assembly_graph)

    def assembly_space(self, **limits) -> AssemblySpace:
        """Create a reusable enumeration for counting, selecting and sampling."""
        return AssemblySpace(self.reaction_graph, self.root_enc, input_heavy_atoms=self.input_heavy_atoms, **limits)

    def count_unique_assemblies(self, *, min_coverage: float | None = 0.0, **limits) -> int:
        """Count all views; pass None for maximum coverage, or a minimum fraction."""
        return self.assembly_space(**limits).count(min_coverage=min_coverage)

    def sample_assemblies(
        self, n: int = 1, *, seed: int | None = None, min_coverage: float | None = None, **limits,
    ) -> list[AssemblyGraph]:
        """Sample maximum-coverage views, or explicitly allow a lower minimum."""
        return self.assembly_space(**limits).sample(n, seed=seed, min_coverage=min_coverage)

    @cached_property
    def input_heavy_atoms(self) -> int:
        """All submitted heavy atoms, including fragments removed in preparation."""
        return smiles_to_mol(self.submission.smiles).GetNumHeavyAtoms()

    def calculate_coverage(self, assembly_graph: AssemblyGraph | None = None) -> float:
        """Identified assembly heavy atoms / ALL submitted heavy atoms.

        Supply an assembly to measure a particular sampled resolution. Without
        one, use the default projection, never the union of alternative routes.
        """
        input_heavy_atoms = self.input_heavy_atoms
        if not input_heavy_atoms:
            return 0.0
        root_heavy_tags = {
            atom.GetIsotope() for atom in self.root.mol.GetAtoms()
            if atom.GetAtomicNum() > 1 and atom.GetIsotope() != 0
        }
        assembly = self.assembly_graph if assembly_graph is None else assembly_graph
        identified_tags: set[int] = set()
        for node in assembly.monomer_nodes():
            if node.is_identified:
                identified_tags.update(get_tags_mol(node.mol))
        return len(identified_tags & root_heavy_tags) / input_heavy_atoms

    def to_dict(self) -> dict[str, Any]:
        """Serialize source data only, using the compact version 2 schema."""
        return {
            "schema_version": 2,
            "root": {
                "node": sorted(self.reaction_graph.nodes).index(self.root_enc),
                "smiles": self.submission.smiles,
                "name": self.submission.name,
                "keep_stereo": self.submission.keep_stereo,
                "neutralize": self.submission.neutralize,
                "canonicalize_tautomer": self.submission.canonicalize_tautomer,
                "match_stereochemistry": self.match_stereochemistry,
            },
            "props": self.props,
            "reaction_rules_hash": self.reaction_rules_hash,
            "matching_rules_hash": self.matching_rules_hash,
            "reaction_graph": self.reaction_graph.to_dict(),
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "Result":
        """Load version 2. Older result formats are intentionally unsupported."""
        if data.get("schema_version") != 2:
            raise ValueError("Unsupported result schema; reparse the input to produce schema_version 2")
        graph = ReactionGraph.from_dict(data["reaction_graph"])
        root = data["root"]
        enc = ReactionGraph._at(list(graph.nodes), root["node"])
        submission = Submission(
            smiles=root["smiles"], name=root["name"], props=data["props"],
            keep_stereo=root["keep_stereo"], neutralize=root["neutralize"],
            canonicalize_tautomer=root["canonicalize_tautomer"],
        )
        # The saved root is authoritative for atom numbering and chemistry.
        mol = Mol(graph.nodes[enc].mol)
        object.__setattr__(submission, "mol", mol)
        object.__setattr__(submission, "inchikey", mol_to_inchikey(remove_tags(mol, in_place=False)))
        object.__setattr__(submission, "stereo_registry", capture_double_bond_stereo(mol))
        return cls(
            submission, graph,
            reaction_rules_hash=data["reaction_rules_hash"],
            matching_rules_hash=data["matching_rules_hash"],
            match_stereochemistry=root["match_stereochemistry"], root_enc=enc,
        )
