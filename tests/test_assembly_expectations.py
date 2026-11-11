"""Independent graph expectations exercise the integration-test contract."""

from dataclasses import replace

import pytest

from retromol.model.result import Result

from .helpers_retromol import ExpectedAssembly, assert_assemblies, parse_compound
from .test_assembly_sampling import layered_result
from .test_result import make_result


FINE_CHAIN = ExpectedAssembly(
    coverage=1.0,
    nodes=(("C", 1), ("C", 1), ("C", 1), ("C", 1)),
    edges=((0, 1, ("C-C:SINGLE",)), (1, 2, ("C-C:SINGLE",)), (2, 3, ("C-C:SINGLE",))),
)
COARSE_CHAIN = ExpectedAssembly(
    coverage=1.0,
    nodes=(("CC", 2), ("CC", 2)),
    edges=((0, 1, ("C-C:SINGLE",)),),
)
WRONG_BRANCH = replace(
    FINE_CHAIN,
    edges=((0, 1, ("C-C:SINGLE",)), (0, 2, ("C-C:SINGLE",)), (0, 3, ("C-C:SINGLE",))),
)


@pytest.mark.parametrize("expected", [
    [FINE_CHAIN, COARSE_CHAIN],
    [COARSE_CHAIN, FINE_CHAIN],
])
def test_all_expected_graphs_match_while_extra_parsed_assemblies_are_allowed(expected):
    result, *_ = layered_result()
    assert result.assembly_space().count() == 4
    assert_assemblies(result, expected)


@pytest.mark.parametrize("expected, missing_index", [
    ([WRONG_BRANCH, FINE_CHAIN], 1),
    ([FINE_CHAIN, WRONG_BRANCH], 2),
])
def test_one_matching_graph_does_not_hide_a_missing_expected_graph(expected, missing_index):
    result, *_ = layered_result()
    with pytest.raises(AssertionError, match="Missing 1 of 2 expected assembly graphs") as exc:
        assert_assemblies(result, expected)
    assert f"#{missing_index}:" in str(exc.value)
    assert f"#{3 - missing_index}:" not in str(exc.value)


def test_nondefault_assembly_matches_without_using_default_selection(monkeypatch):
    result, *_ = layered_result()
    assert len(result.assembly_graph.monomer_nodes()) == 4

    def unexpected_default(_):
        raise AssertionError("Expected graphs must search all valid assemblies")

    monkeypatch.setattr(Result, "assembly_graph", property(unexpected_default))
    assert_assemblies(result, [COARSE_CHAIN])


def test_fixture_node_indices_and_edge_orientation_do_not_matter():
    result, *_ = layered_result()
    # A different numbering of the same four-node path.
    shuffled = replace(
        FINE_CHAIN,
        edges=((2, 0, ("C-C:SINGLE",)), (3, 2, ("C-C:SINGLE",)), (1, 3, ("C-C:SINGLE",))),
    )
    assert_assemblies(result, [shuffled])


@pytest.mark.parametrize("expected", [
    WRONG_BRANCH,
    replace(FINE_CHAIN, edges=((0, 1, ("C-C:DOUBLE",)), *FINE_CHAIN.edges[1:])),
    replace(FINE_CHAIN, edges=FINE_CHAIN.edges[:-1]),
    replace(FINE_CHAIN, nodes=(("N", 1), *FINE_CHAIN.nodes[1:])),
    replace(FINE_CHAIN, nodes=(("C", 2), *FINE_CHAIN.nodes[1:])),
    replace(FINE_CHAIN, nodes=FINE_CHAIN.nodes[:3], edges=FINE_CHAIN.edges[:2]),
    replace(FINE_CHAIN, coverage=0.75),
], ids=["wrong-topology", "wrong-bond-order", "missing-edge", "wrong-label",
        "wrong-atom-count", "missing-repeated-monomer", "wrong-coverage"])
def test_graph_mismatches_fail_even_when_monomer_names_overlap(expected):
    result, *_ = layered_result()
    with pytest.raises(AssertionError, match="Missing 1 of 1 expected assembly graphs"):
        assert_assemblies(result, [expected])


def test_unresolved_regions_and_their_connections_are_part_of_the_graph():
    result = make_result("CCO", [("[1CH4]", True), ("[2CH3][3OH]", False)])
    expected = ExpectedAssembly(
        coverage=1 / 3,
        nodes=(("C", 1), (None, 2)),
        edges=((0, 1, ("C-C:SINGLE",)),),
    )
    assert_assemblies(result, [expected])
    with pytest.raises(AssertionError, match="Missing 1 of 1 expected assembly graphs"):
        assert_assemblies(result, [replace(expected, nodes=(("C", 1),), edges=())])


def test_cyclic_dipeptide_requires_both_root_bonds(ruleset):
    result = parse_compound("CN1C(C)C(=O)N(C)C(C)C1=O", ruleset)
    # Two residues, six original heavy atoms each, joined by TWO peptide bonds.
    expected = ExpectedAssembly(
        coverage=1.0,
        nodes=(("N-methylalanine", 6), ("N-methylalanine", 6)),
        edges=((0, 1, ("C-N:SINGLE", "C-N:SINGLE")),),
    )
    assert_assemblies(result, [expected])
    with pytest.raises(AssertionError, match="Missing 1 of 1 expected assembly graphs"):
        assert_assemblies(result, [replace(expected, edges=((0, 1, ("C-N:SINGLE",)),))])


def test_norleucine_peptide_has_two_residues_and_one_peptide_bond(ruleset):
    result = parse_compound("CC(N)C(=O)NC(CCCC)C(=O)O", ruleset)
    expected = ExpectedAssembly(
        coverage=1.0,
        nodes=(("alanine", 5), ("norleucine", 9)),
        edges=((0, 1, ("C-N:SINGLE",)),),
    )
    assert_assemblies(result, [expected])
    assert_assemblies(Result.from_dict(result.to_dict()), [expected])


def test_stereo_qualified_node_labels_are_not_interchangeable(stereo_ruleset):
    result = parse_compound("N[C@@H](C)C(=O)O", stereo_ruleset)
    expected = ExpectedAssembly(coverage=1.0, nodes=(("alanine^L", 6),), edges=())
    assert_assemblies(result, [expected])
    with pytest.raises(AssertionError, match="Missing 1 of 1 expected assembly graphs"):
        assert_assemblies(result, [replace(expected, nodes=(("alanine^D", 6),))])


def test_empty_expected_graph_list_is_an_error():
    result, *_ = layered_result()
    with pytest.raises(ValueError, match="At least one expected assembly"):
        assert_assemblies(result, [])
