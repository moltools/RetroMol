"""Assembly alternatives, compact persistence, pruning and provenance."""

import json
from dataclasses import replace

import pytest
import yaml

from retromol.chem.mol import smiles_to_mol
from retromol.chem.tagging import get_tags_mol
from retromol.model.assembly_sampling import AssemblyConstraintError, AssemblyLimitError, AssemblySpace, count_unique_assemblies, sample_assemblies
from retromol.model.reaction_graph import ReactionGraph, ReactionStep
from retromol.model.readout import LinearReadout
from retromol.model.result import Result
from retromol.model.rules import MatchingRule, ReactionRule, RuleSet
from retromol.model.submission import Submission
from retromol.pipelines import parsing


STEP = ReactionStep('contested', ('split',), ('rule-1',))


def add(graph, smiles, *, known=False, terminal=False):
    enc = graph.add_node(smiles_to_mol(smiles))
    if known:
        graph.nodes[enc].identify([MatchingRule(
            name=graph.nodes[enc].smiles, smiles=graph.nodes[enc].smiles,
            props={}, pseudonyms=[], terminal=terminal,
        )])
    else:
        graph.nodes[enc].identify([])
    return enc


def edge(graph, source, *children, step=STEP):
    graph.add_edge(source, [graph.nodes[c].mol for c in children], step)


def layered_result():
    submission = Submission('CCCC')
    graph = ReactionGraph()
    root = graph.add_node(submission.mol)
    left = add(graph, '[1CH3][2CH3]', known=True)
    right = add(graph, '[3CH3][4CH3]', known=True)
    singles = [add(graph, f'[{tag}CH4]', known=True, terminal=True) for tag in range(1, 5)]
    edge(graph, root, left, right)
    edge(graph, left, *singles[:2])
    edge(graph, right, *singles[2:])
    return Result(submission, graph), left, right, singles


def signature(assembly):
    nodes = tuple(sorted(assembly.monomer_ids()))
    bonds = sorted(tuple(sorted((link.a1_tag, link.a2_tag)))
                   for _, _, links in assembly.edges_with_bonds() for link in links)
    return nodes, bonds


def test_independent_regions_allow_mixed_resolution_and_rebuild_root_bonds():
    result, left, right, singles = layered_result()
    space = result.assembly_space()
    assert space.count() == 4
    assert set(space.frontiers) == {
        frozenset([left, right]), frozenset([*singles[:2], right]),
        frozenset([left, *singles[2:]]), frozenset(singles),
    }
    coarse = space.build(frozenset([left, right]))
    assert signature(coarse)[1] == [(2, 3)]
    fine = space.build(frozenset(singles))
    assert signature(fine)[1] == [(1, 2), (2, 3), (3, 4)]
    for assembly in space.sample(10, seed=17):
        assert result.calculate_coverage(assembly) == 1
        assert assembly.g.nodes['unassigned']['tags'] == set()
        tags = [get_tags_mol(n.mol) for n in assembly.monomer_nodes()]
        assert sum(map(len, tags)) == len(set.union(*tags)) == 4


def test_duplicate_routes_and_product_permutations_count_once():
    result, left, right, singles = layered_result()
    edge(result.reaction_graph, result.root_enc, right, left,
         step=ReactionStep('contested', ('other route',), ('rule-2',)))
    edge(result.reaction_graph, result.root_enc, *singles)
    assert result.count_unique_assemblies() == 4
    assert len({signature(a)[0] for a in result.sample_assemblies(100, seed=42)}) == 4


def terminal_bypass_result():
    result, left, right, singles = layered_result()
    graph = result.reaction_graph
    node = graph.nodes[left]
    graph.nodes[left] = replace(node, identity=replace(
        node.identity, matched_rule=replace(node.identity.matched_rule, terminal=True),
    ))
    # This shortcut decomposes the terminal left unit before recognizing it.
    edge(graph, result.root_enc, *singles)
    return result, left, right, singles


def test_terminal_boundaries_apply_across_routes_before_coverage_filtering():
    result, left, right, singles = terminal_bypass_result()
    payload = json.loads(json.dumps(result.to_dict()))
    space = result.assembly_space()
    assert set(space.frontiers) == {frozenset([left, right]), frozenset([left, *singles[2:]])}
    assert space.count() == space.count(min_coverage=None) == 2
    for minimum in [None, 0, .5, 1]:
        samples = space.sample(10, seed=42, min_coverage=minimum)
        assert len(samples) == 2
        assert all(left in assembly.monomer_ids() for assembly in samples)
        assert all(result.calculate_coverage(assembly) == 1 for assembly in samples)
    with pytest.raises(ValueError, match='not a valid assembly frontier'):
        space.build(frozenset(singles))
    loaded = Result.from_dict(payload)
    assert [signature(a) for a in loaded.sample_assemblies(10, seed=42)] == [
        signature(a) for a in result.sample_assemblies(10, seed=42)
    ]
    assert json.loads(json.dumps(result.to_dict())) == payload


def test_default_assembly_and_readout_replace_a_cheaper_terminal_bypass():
    result, left, right, singles = terminal_bypass_result()
    graph = result.reaction_graph
    old_default = parsing.extract_min_edge_synthesis_subgraph(
        graph, result.root_enc, edge_base_cost=.25, nonterminal_leaf_penalty=100,
    ).graph
    assert {n.enc for n in old_default.get_leaf_nodes()} == set(singles)
    expected = {left, *singles[2:]}
    assert set(result.assembly_graph.monomer_ids()) == expected
    assert set(result.linear_readout.assembly_graph.monomer_ids()) == expected
    assert {n.enc for n in result.selected_reaction_graph.get_leaf_nodes()} == expected
    assert result.calculate_coverage() == 1
    assert all(e in graph.edges for e in result.selected_reaction_graph.edges)


def test_nonterminal_setting_explicitly_allows_finer_views_again():
    result, left, _, singles = terminal_bypass_result()
    graph = result.reaction_graph
    node = graph.nodes[left]
    graph.nodes[left] = replace(node, identity=replace(
        node.identity, matched_rule=replace(node.identity.matched_rule, terminal=False),
    ))
    assert result.count_unique_assemblies() == 4
    assert frozenset(singles) in result.assembly_space().frontiers


def test_crossing_terminal_assignments_compete_but_subdivisions_are_excluded():
    submission = Submission('CCC')
    graph = ReactionGraph()
    root = graph.add_node(submission.mol)
    ab = add(graph, '[1CH3][2CH3]', known=True, terminal=True)
    bc = add(graph, '[2CH3][3CH3]', known=True, terminal=True)
    a, b, c = [add(graph, f'[{tag}CH4]', known=True, terminal=True) for tag in range(1, 4)]
    edge(graph, root, ab, c)
    edge(graph, root, a, bc)
    edge(graph, root, a, b, c)
    result = Result(submission, graph)
    space = result.assembly_space()
    assert set(space.frontiers) == {frozenset([ab, c]), frozenset([a, bc])}
    assert space.count() == 2
    assert not space.preserves_terminals(frozenset([a, b, c]))
    assert space.preserves_terminals(frozenset(result.assembly_graph.monomer_ids()))


def test_disconnected_terminal_match_does_not_constrain_another_root():
    result, *_ = layered_result()
    add(result.reaction_graph, '[1CH3][2CH2][3CH3]', known=True, terminal=True)
    assert result.count_unique_assemblies() == 4


def test_terminal_unit_cannot_be_hidden_by_losing_some_of_its_original_atoms():
    result, left, right, singles = terminal_bypass_result()
    # Dropping the second carbon is not an acceptable way to avoid the boundary.
    edge(result.reaction_graph, result.root_enc, singles[0], right)
    assert result.count_unique_assemblies() == 2
    assert all(left in a.monomer_ids() for a in result.sample_assemblies(10, min_coverage=0))


def test_terminal_boundaries_ignore_reaction_added_atoms():
    submission = Submission('CCC')
    graph = ReactionGraph()
    root = graph.add_node(submission.mol)
    known = add(graph, '[1CH3][2CH2][99OH]', known=True, terminal=True)
    singles = [add(graph, f'[{tag}CH4]', known=True, terminal=True) for tag in range(1, 4)]
    edge(graph, root, known, singles[2])
    edge(graph, root, *singles)
    result = Result(submission, graph)
    assert result.count_unique_assemblies() == 1
    assert set(result.sample_assemblies()[0].monomer_ids()) == {known, singles[2]}


def test_incompatible_terminal_constraints_do_not_invent_a_route_or_relax_protection():
    result, left, right, singles = layered_result()
    graph = result.reaction_graph
    for enc in [left, right]:
        node = graph.nodes[enc]
        graph.nodes[enc] = replace(node, identity=replace(
            node.identity, matched_rule=replace(node.identity.matched_rule, terminal=True),
        ))
    # Each route preserves one terminal unit and splits the other. The graph
    # never explored a route that keeps both intact, so we must not invent one.
    graph.edges.clear()
    graph.out_edges = {enc: [] for enc in graph.nodes}
    edge(graph, result.root_enc, left, *singles[2:])
    edge(graph, result.root_enc, *singles[:2], right)
    assert result.count_unique_assemblies() == 0
    assert result.sample_assemblies(10, min_coverage=0) == []
    with pytest.raises(AssemblyConstraintError, match='preserves the identified terminal units'):
        _ = result.assembly_graph


def test_equal_names_at_different_atom_positions_remain_distinct():
    submission = Submission('CCC')
    graph = ReactionGraph()
    root = graph.add_node(submission.mol)
    a = add(graph, '[1CH4]', known=True)
    bc = add(graph, '[2CH3][3CH3]', known=True)
    ab = add(graph, '[1CH3][2CH3]', known=True)
    c = add(graph, '[3CH4]', known=True)
    edge(graph, root, a, bc)
    edge(graph, root, ab, c)
    space = AssemblySpace(graph, root)
    assert space.count() == 2
    assert len({signature(x)[0] for x in space.sample(10, seed=1)}) == 2


def test_overlapping_products_are_not_compatible_assemblies():
    result, left, right, singles = layered_result()
    edge(result.reaction_graph, result.root_enc, left, singles[0], right)
    assert result.count_unique_assemblies() == 4
    with pytest.raises(ValueError, match='not a valid'):
        result.assembly_space().build(frozenset([left, singles[0], right]))


def test_reaction_introduced_byproducts_do_not_create_extra_assemblies():
    result, left, right, _ = layered_result()
    water = add(result.reaction_graph, 'O', known=True)
    edge(result.reaction_graph, result.root_enc, left, right, water)
    assert result.count_unique_assemblies() == 4
    assert all(water not in f for f in result.assembly_space().frontiers)


def test_seeded_sampling_is_reproducible_after_json_roundtrip():
    result, *_ = layered_result()
    loaded = Result.from_dict(json.loads(json.dumps(result.to_dict())))
    expected = [signature(a) for a in result.sample_assemblies(3, seed=13)]
    assert [signature(a) for a in loaded.sample_assemblies(3, seed=13)] == expected
    assert count_unique_assemblies(loaded.reaction_graph, loaded.root_enc) == 4
    assert len(sample_assemblies(loaded.reaction_graph, loaded.root_enc, 10, seed=13)) == 4
    assert result.sample_assemblies(0) == []
    with pytest.raises(ValueError, match='nonnegative integer'):
        result.sample_assemblies(-1)


def partially_identified_result():
    result, _, _, singles = layered_result()
    graph = result.reaction_graph
    graph.nodes[singles[0]] = replace(graph.nodes[singles[0]], identified=False, identity=None)
    return result


def test_default_sampling_uses_global_maximum_not_default_projection_coverage():
    result = partially_identified_result()
    # The cost-based fine-grained projection is not necessarily the best-covered.
    assert result.calculate_coverage() == .75
    space = result.assembly_space()
    assert space.count() == 4
    assert space.max_coverage == 1
    assert space.count(min_coverage=None) == 2
    assert result.count_unique_assemblies(min_coverage=None) == 2
    best = space.sample(10, seed=42)
    assert len(best) == len({signature(a)[0] for a in best}) == 2
    assert all(result.calculate_coverage(a) == 1 for a in best)
    for minimum, expected_count in [(0, 4), (.75, 4), (.76, 2), (1, 2)]:
        samples = space.sample(10, seed=42, min_coverage=minimum)
        assert len(samples) == space.count(min_coverage=minimum) == expected_count
        assert all(result.calculate_coverage(a) >= minimum for a in samples)


def test_maximum_coverage_below_one_includes_all_input_atoms_and_roundtrips():
    result = partially_identified_result()
    result = replace(result, submission=replace(result.submission, smiles='CCCC.[Na+]'))
    payload = json.loads(json.dumps(result.to_dict()))
    loaded = Result.from_dict(payload)
    space = loaded.assembly_space()
    assert space.max_coverage == .8  # largest fragment has 4 of the 5 input heavy atoms
    assert len(space.sample(10, seed=1)) == 2
    assert all(loaded.calculate_coverage(a) == .8 for a in space.sample(10, seed=1))
    assert space.count(min_coverage=.81) == 0
    assert space.sample(10, min_coverage=.81) == []  # never backfill with lower coverage
    assert [signature(a) for a in loaded.sample_assemblies(1, seed=17)] == [
        signature(a) for a in result.sample_assemblies(1, seed=17)
    ]
    assert json.loads(json.dumps(loaded.to_dict())) == payload
    assert count_unique_assemblies(loaded.reaction_graph, loaded.root_enc,
                                   input_heavy_atoms=5, min_coverage=.81) == 0
    assert sample_assemblies(loaded.reaction_graph, loaded.root_enc, 10,
                             input_heavy_atoms=5, min_coverage=.81) == []


@pytest.mark.parametrize('minimum', [-.1, 1.1, float('nan'), float('inf'), -float('inf')])
def test_sampling_rejects_invalid_coverage_thresholds(minimum):
    result, *_ = layered_result()
    with pytest.raises(ValueError, match='finite fraction between 0 and 1'):
        result.sample_assemblies(10, min_coverage=minimum)
    with pytest.raises(ValueError, match='finite fraction between 0 and 1'):
        result.count_unique_assemblies(min_coverage=minimum)


def test_zero_coverage_still_has_a_best_assembly():
    submission = Submission('CCO')
    graph = ReactionGraph()
    graph.add_node(submission.mol)
    result = Result(submission, graph)
    space = result.assembly_space()
    assert space.max_coverage == 0
    assert space.count(min_coverage=None) == 1
    assert result.calculate_coverage(space.sample()[0]) == 0
    assert space.sample(min_coverage=.01) == []


def test_count_and_sample_raise_instead_of_silently_truncating():
    result, *_ = layered_result()
    with pytest.raises(AssemblyLimitError, match='max_assemblies=2'):
        result.count_unique_assemblies(max_assemblies=2)
    with pytest.raises(AssemblyLimitError, match='max_assemblies=2'):
        result.sample_assemblies(1, seed=1, max_assemblies=2)
    with pytest.raises(AssemblyLimitError, match='max_combinations=1'):
        result.count_unique_assemblies(max_combinations=1)
    with pytest.raises(ValueError, match='positive'):
        result.assembly_space(max_assemblies=0)


def test_pruning_preserves_unknown_boundary_and_identified_siblings():
    submission = Submission('CCCO')
    graph = ReactionGraph()
    root = graph.add_node(submission.mol)
    known = add(graph, '[1CH4]', known=True)
    unknown = add(graph, '[2CH3][3CH2][4OH]')
    tail = add(graph, '[2CH3][3CH]=[4O]')
    leaves = [add(graph, '[2CH4]'), add(graph, '[3CH2]=[4O]')]
    edge(graph, root, known, unknown)
    edge(graph, unknown, tail)
    edge(graph, tail, *leaves)
    pruned = graph.prune_unidentified(root)
    assert set(pruned.nodes) == {root, known, unknown}
    assert pruned.nodes[unknown].pruned
    assert not graph.nodes[unknown].pruned  # source graph was not mutated
    result = Result(submission, pruned)
    assert result.count_unique_assemblies() == 1
    assembly = result.sample_assemblies(seed=3)[0]
    assert set.union(*(get_tags_mol(n.mol) for n in assembly.monomer_nodes())) == {1, 2, 3, 4}
    assert result.calculate_coverage(assembly) == .25
    restored = Result.from_dict(result.to_dict())
    assert restored.reaction_graph.nodes[unknown].pruned


def test_wholly_unidentified_graph_collapses_to_root():
    submission = Submission('CCO')
    graph = ReactionGraph()
    root = graph.add_node(submission.mol)
    child = add(graph, '[1CH3][2CH]=[3O]')
    edge(graph, root, child)
    edge(graph, child, root)
    result = Result(submission, graph.prune_unidentified(root))
    assert set(result.reaction_graph.nodes) == {root}
    assert result.reaction_graph.edges == []
    assert result.root.pruned
    assert result.count_unique_assemblies() == 1
    assert result.calculate_coverage() == 0


def test_pruning_drops_wholly_unproductive_alternative_instead_of_making_it_cheaper():
    result, *_ = layered_result()
    graph = result.reaction_graph
    unknown = add(graph, '[1CH2]=[2CH][3CH2][4CH3]')
    deeper = add(graph, '[1CH2]=[2CH][3CH]=[4CH2]')
    edge(graph, result.root_enc, unknown)
    edge(graph, unknown, deeper)
    pruned = graph.prune_unidentified(result.root_enc)
    assert unknown not in pruned.nodes and deeper not in pruned.nodes
    assert len(pruned.out_edges[result.root_enc]) == 1
    assert Result(result.submission, pruned).count_unique_assemblies() == 4


def test_identified_parent_drops_wholly_unknown_descendants():
    graph = ReactionGraph()
    root = add(graph, '[1CH4]', known=True)
    child = add(graph, '[1CH3]O')
    edge(graph, root, child)
    pruned = graph.prune_unidentified(root)
    assert set(pruned.nodes) == {root}
    assert pruned.nodes[root].is_identified and pruned.nodes[root].pruned
    assert AssemblySpace(pruned, root).frontiers == (frozenset([root]),)


def test_collapsed_unknown_region_does_not_become_a_cheap_default_shortcut():
    result, _, _, singles = layered_result()
    graph = result.reaction_graph
    unknown = add(graph, '[2CH3][3CH2][4CH3]')
    graph.nodes[unknown] = replace(graph.nodes[unknown], pruned=True)
    edge(graph, result.root_enc, singles[0], unknown)
    costs = dict(edge_base_cost=4, nonterminal_leaf_penalty=100)
    plain = parsing.extract_min_edge_synthesis_subgraph(graph, result.root_enc, **costs)
    weighted = parsing.extract_min_edge_synthesis_subgraph(
        graph, result.root_enc, weight_pruned_atoms=True, **costs,
    )
    assert not plain.solved
    assert weighted.solved
    assert {n.enc for n in weighted.graph.get_leaf_nodes()} == set(singles)


def test_productive_cycles_terminate_and_keep_identified_exit():
    graph = ReactionGraph()
    root = add(graph, '[1CH4]')
    intermediate = add(graph, '[1CH3]O')
    identified = add(graph, '[1CH3]N', known=True, terminal=True)
    edge(graph, root, intermediate)
    edge(graph, intermediate, root)
    edge(graph, intermediate, identified)
    pruned = graph.prune_unidentified(root)
    assert len(pruned.edges) == 3
    assert AssemblySpace(pruned, root).frontiers == (frozenset([identified]),)
    edge(graph, root, root, identified)
    assert len(graph.prune_unidentified(root).edges) == 3


def test_compact_schema_interns_records_and_keeps_input_options_and_props():
    result, *_ = layered_result()
    result = replace(result, submission=Submission(
        'CCCC.[Na+]', name='example', props={'source': 1, 'root': 'user property'},
        keep_stereo=False, neutralize=False, canonicalize_tautomer=True,
    ), reaction_rules_hash='a' * 64, matching_rules_hash='b' * 64, match_stereochemistry=True)
    # Derive views first; these caches must still never enter the saved result.
    _ = result.linear_readout
    data = json.loads(json.dumps(result.to_dict()))
    assert set(data) == {'schema_version', 'root', 'props', 'reaction_graph',
                         'reaction_rules_hash', 'matching_rules_hash'}
    g = data['reaction_graph']
    assert len(g['identities']) == 2  # ethane and methane, each shared across occurrences
    assert len(g['steps']) == 1
    assert 'out_edges' not in g
    assert all('enc' not in n and 'tagged_smiles' not in n for n in g['nodes'])
    loaded = Result.from_dict(data)
    assert loaded.root_enc == result.root_enc
    assert loaded.props == {'source': 1, 'root': 'user property'}
    assert loaded.submission.name == 'example'
    assert not loaded.submission.keep_stereo and not loaded.submission.neutralize
    assert loaded.submission.canonicalize_tautomer and loaded.match_stereochemistry
    assert loaded.reaction_rules_hash == 'a' * 64
    assert loaded.matching_rules_hash == 'b' * 64
    assert loaded.calculate_coverage() == .8
    assert loaded.count_unique_assemblies() == 4
    assert 'linear_readout' not in loaded.__dict__
    assert json.loads(json.dumps(loaded.to_dict())) == data


@pytest.mark.parametrize('index', [-1, 1_000_000, True])
def test_invalid_serialized_node_references_are_rejected(index):
    result, *_ = layered_result()
    data = result.to_dict()
    data['root']['node'] = index
    with pytest.raises(ValueError, match='Invalid graph reference'):
        Result.from_dict(data)


def test_old_result_format_is_explicitly_rejected():
    with pytest.raises(ValueError, match='reparse'):
        Result.from_dict({'submission': {}, 'reaction_graph': {}, 'linear_readout': {}})


def test_parse_keeps_competing_routes_and_does_not_compute_readouts(monkeypatch):
    result, left, right, singles = layered_result()
    edge(result.reaction_graph, result.root_enc, *singles)
    rules = RuleSet(False, [], [])
    monkeypatch.setattr(parsing, 'process_mol', lambda *args: result.reaction_graph)

    def unexpected_selection(*args, **kwargs):
        raise AssertionError('Parsing and saving must not choose an assembly')

    monkeypatch.setattr(parsing, 'extract_min_edge_synthesis_subgraph', unexpected_selection)
    monkeypatch.setattr(LinearReadout, 'from_assembly_graph', unexpected_selection)
    parsed = parsing.run_retromol(result.submission, rules)
    assert len(parsed.reaction_graph.out_edges[parsed.root_enc]) == 2
    assert parsed.count_unique_assemblies() == 4
    assert parsed.to_dict()['reaction_rules_hash'] == rules.reaction_rules_hash
    assert parsed.matching_rules_hash == rules.matching_rules_hash


def test_decanoic_lipopeptide_has_five_resolutions(ruleset):
    result = parsing.run_retromol(Submission('CCCCCCCCCC(=O)NCC(=O)NC(C)C(=O)O'), ruleset)
    space = result.assembly_space()
    assert space.count() == 5
    assemblies = space.sample(10, seed=42)
    assert sorted(len(a.monomer_ids()) for a in assemblies) == [3, 4, 5, 6, 7]
    fatty_names = {n.identity.name for a in assemblies for n in a.monomer_nodes()
                   if 'acid' in n.identity.name}
    assert fatty_names == {'decanoic acid', 'octanoic acid', 'hexanoic acid', 'butanoic acid', 'acetic acid'}
    assert all(result.calculate_coverage(a) == 1 for a in assemblies)
    sequences = [LinearReadout.from_assembly_graph(a).primary_sequence() for a in assemblies]
    assert sorted(map(len, sequences)) == [3, 4, 5, 6, 7]


def test_rule_hashes_include_definition_metadata_and_order(tmp_path):
    reaction = ReactionRule('cut', '[C:1]-[C:2]>>[C:1].[C:2]', {'a': 1, 'b': 2})
    matching = MatchingRule('methane', 'C', {}, [], terminal=False)
    rules = RuleSet(False, [reaction], [matching])
    changed_rxn = RuleSet(False, [replace(reaction, allowed_in_bulk=True)], [matching])
    changed_match = RuleSet(False, [reaction], [replace(matching, terminal=True)])
    assert rules.reaction_rules_hash != changed_rxn.reaction_rules_hash
    assert rules.matching_rules_hash == changed_rxn.matching_rules_hash
    assert rules.matching_rules_hash != changed_match.matching_rules_hash
    assert rules.reaction_rules_hash == changed_match.reaction_rules_hash
    assert rules.reaction_rules_hash == RuleSet(False, [replace(reaction, props={'b': 2, 'a': 1})], []).reaction_rules_hash
    a = replace(matching, name='first')
    b = replace(matching, name='second')
    assert RuleSet(False, [], [a, b]).matching_rules_hash != RuleSet(False, [], [b, a]).matching_rules_hash
    rxn_path, mxn_path = tmp_path / 'rxn.yml', tmp_path / 'mxn.yml'
    rxn_path.write_text(yaml.safe_dump([reaction.to_dict()]))
    mxn_path.write_text(yaml.safe_dump([matching.to_dict()]))
    first = RuleSet.load_from_files(rxn_path, mxn_path)
    rxn_path.write_text('# different formatting\n' + yaml.safe_dump([reaction.to_dict()], default_flow_style=True))
    second = RuleSet.load_from_files(rxn_path, mxn_path)
    assert first.reaction_rules_hash == second.reaction_rules_hash == rules.reaction_rules_hash
    assert first.matching_rules_hash == second.matching_rules_hash == rules.matching_rules_hash
