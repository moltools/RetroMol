"""Regression tests for coverage of identified assembly monomers."""

from dataclasses import replace

import pytest

from retromol.chem.mol import mol_to_smiles, smiles_to_mol
from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.reaction_graph import ReactionGraph, ReactionStep
from retromol.model.readout import LinearReadout
from retromol.model.result import Result
from retromol.model.rules import MatchingRule
from retromol.model.submission import Submission


def make_result(smiles: str, fragments: list[tuple[str, bool]]) -> Result:
    """Build an identified parent with explicitly identified/unknown assembly leaves."""
    submission = Submission(smiles)
    graph = ReactionGraph()
    root = graph.add_node(submission.mol)
    parent_rule = MatchingRule(
        name="identified parent", smiles=mol_to_smiles(submission.mol, include_tags=True),
        props={}, pseudonyms=[], terminal=False,
    )
    graph.nodes[root].identify([parent_rule])
    children = graph.add_edge(
        root,
        [smiles_to_mol(smi) for smi, _ in fragments],
        ReactionStep(kind="contested", names=("split",), rule_ids=("split",)),
    )
    for enc, (tagged_smiles, identified) in zip(children, fragments):
        node = graph.nodes[enc]
        rules = [MatchingRule(name=node.smiles, smiles=tagged_smiles, props={}, pseudonyms=[])] if identified else []
        node.identify(rules)
    return Result(submission, graph)


def test_coverage_ignores_identified_intermediates_and_unknown_leaves() -> None:
    result = make_result("CCOC", [("[1CH3][2CH3]", True), ("[3OH][4CH3]", False)])

    assert len(result.reaction_graph.identified_nodes) == 2
    assert result.calculate_coverage() == pytest.approx(0.5)
    assert Result.from_dict(result.to_dict()).calculate_coverage() == pytest.approx(0.5)


def test_coverage_counts_assembly_monomers_without_linear_paths(monkeypatch) -> None:
    result = make_result("CCOC", [("[1CH3][2CH3]", True), ("[3OH][4CH3]", True)])

    def unexpected_readout(*args, **kwargs):
        raise AssertionError("Coverage must not request sequences")

    monkeypatch.setattr(LinearReadout, "from_assembly_graph", unexpected_readout)
    assert result.calculate_coverage() == 1.0


def test_coverage_keeps_atoms_missing_from_assembly_in_denominator() -> None:
    result = make_result("CCOC", [("[1CH3][2CH3]", True)])

    assert result.calculate_coverage() == pytest.approx(0.5)


def test_coverage_includes_all_submitted_fragments_in_denominator() -> None:
    result = make_result("CCO.[Na+]", [("[1CH3][2CH2][3OH]", True)])

    assert result.submission.mol.GetNumHeavyAtoms() == 3
    assert result.calculate_coverage() == pytest.approx(3 / 4)
    assert Result.from_dict(result.to_dict()).calculate_coverage() == pytest.approx(3 / 4)
    assert result.assembly_space().max_coverage == pytest.approx(3 / 4)
    assert result.sample_assemblies(min_coverage=1) == []


def test_coverage_ignores_hydrogens_and_reaction_introduced_atoms() -> None:
    # The isotope-labelled H survives SMILES parsing and gets root tag 1.
    # Carbon tags 2/3 are covered; oxygen tag 4 is absent. Tag 99 and the
    # untagged oxygen were introduced by a reaction and cannot count as input.
    result = make_result("[2H]CCO", [("[1H][2CH2][3CH2][99CH2]O", True)])

    assert result.submission.mol.GetNumAtoms() == 4
    assert result.calculate_coverage() == pytest.approx(2 / 3)
    space = result.assembly_space()
    assert {result.calculate_coverage(a) for a in space.sample(10, min_coverage=0)} == {2 / 3, 1}
    assert space.count(min_coverage=.7) == 1


def test_coverage_does_not_drop_untagged_input_atoms_from_denominator() -> None:
    result = make_result("CCO", [("[1CH3][2CH2]O", True)])
    result.root.mol.GetAtomWithIdx(2).SetIsotope(0)

    assert result.calculate_coverage() == pytest.approx(2 / 3)


@pytest.mark.parametrize("smiles", ["[2H][2H]", ""])
def test_coverage_with_no_heavy_atoms_is_zero(smiles: str) -> None:
    result = make_result(smiles, [])

    assert result.calculate_coverage() == 0.0
    assert result.assembly_space().max_coverage == 0.0
    assert result.calculate_coverage(result.sample_assemblies()[0]) == 0.0


def test_coverage_of_empty_assembly_is_zero() -> None:
    result = make_result("CCO", [])

    assert result.calculate_coverage(AssemblyGraph.build(result.root.mol, [])) == 0.0


def test_coverage_counts_repeated_atom_tags_only_once() -> None:
    result = make_result("CCO", [("[1CH3][2CH3]", True)])
    graph = result.assembly_graph.g
    node = result.assembly_graph.monomer_nodes()[0]
    graph.add_node(
        "overlapping monomer", tags=graph.nodes[node.enc]["tags"],
        identity=node.identity, molnode=replace(node, enc="overlapping monomer"),
    )

    assert result.calculate_coverage() == pytest.approx(2 / 3)
