"""Linear assembly parts, role boundaries, and preservation of branched units."""

import json
import random
from dataclasses import replace

import networkx as nx
import pytest

from retromol.chem.mol import smiles_to_mol
from retromol.model.assembly_graph import AssemblyGraph, RootBondLink
from retromol.model.identity import MolIdentity
from retromol.model.reaction_graph import MolNode
from retromol.model.readout import LinearReadout, assembly_to_primary_sequence, split_named_paths
from retromol.model.result import Result
from retromol.model.rules import MatchingRule, RuleSet
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol


def assembly(names, edges=(), roles=None):
    """Synthetic connectivity fixture; graph IDs deliberately differ from enc."""
    g = nx.Graph()
    for i, name in enumerate(names):
        role = roles[i] if roles else None
        identity = None if name is None else MolIdentity(MatchingRule(
            name, 'C', {}, [], sequence_role=role,
        ))
        node = MolNode(f'mol-{i}', smiles_to_mol(f'[{i + 1}CH4]'), 'C', identity, name is not None)
        g.add_node(i, molnode=node, tags={i + 1}, identity=identity)
    for u, v, elements in edges:
        a, b = elements.split('-')
        bond = RootBondLink(u, v, u + 1, v + 1, a, b, 'SINGLE', 1.0)
        g.add_edge(u, v, bonds=[bond], n_bonds=1)
    return AssemblyGraph(g, validate_upon_initialization=True)


def test_unannotated_rules_keep_the_bond_fallback():
    a = assembly(['A', 'B', 'C', 'sugar'], [
        (0, 1, 'C-N'), (1, 2, 'C-C'), (1, 3, 'C-O'),
    ])
    assert assembly_to_primary_sequence(a) == [['A', 'B', 'C'], ['sugar']]


def test_roles_separate_c_and_n_linked_tailoring_and_group_connected_sugars():
    a = assembly(['A', 'B', 'C', 'methylation', 'sugar1', 'sugar2'], [
        (0, 1, 'C-N'), (1, 2, 'C-C'), (1, 3, 'C-N'),
        (2, 4, 'C-C'), (4, 5, 'C-O'),
    ], roles=[None, 'backbone', None, 'tailoring', 'tailoring', 'tailoring'])
    before = a.to_dict()
    expected = [['A', 'B', 'C'], ['sugar1', 'sugar2'], ['methylation']]
    readout = LinearReadout.from_assembly_graph(a)
    assert readout.primary_sequence_parts() == expected
    assert assembly_to_primary_sequence(a) == expected
    assert split_named_paths(readout.primary_sequence(), key=lambda x: x) == expected
    assert a.to_dict() == before
    restored = AssemblyGraph.from_dict(json.loads(json.dumps(before)))
    assert assembly_to_primary_sequence(restored) == expected
    assert LinearReadout.from_dict(readout.to_dict()).primary_sequence_parts() == expected


@pytest.mark.parametrize('roles, expected', [
    (['backbone', 'backbone'], [['A', 'B']]),
    (['backbone', None], [['A'], ['B']]),
    (['backbone', 'tailoring'], [['A'], ['B']]),
    (['tailoring', 'tailoring'], [['A', 'B']]),
])
def test_explicit_equal_roles_retain_oxygen_links(roles, expected):
    a = assembly(['A', 'B'], [(0, 1, 'C-O')], roles)
    assert assembly_to_primary_sequence(a) == expected


def test_branches_are_retained_once_without_repeating_the_branch_point():
    a = assembly(['A', 'B', 'C', 'D'], [(0, i, 'C-C') for i in (1, 2, 3)])
    assert assembly_to_primary_sequence(a) == [['B', 'A', 'C'], ['D']]


def test_longest_path_compares_all_disconnected_components():
    a = assembly(['A', 'B', 'C', 'D', 'E'], [(0, 1, 'C-C'), (1, 2, 'C-C'), (3, 4, 'C-C')])
    assert [node.identity.name for node in a.longest_path()] == ['A', 'B', 'C']
    assert assembly_to_primary_sequence(a) == [['A', 'B', 'C'], ['D', 'E']]


def test_ring_has_a_stable_cut_and_no_repeated_start():
    a = assembly(['A', 'B', 'C', 'D'], [(i, (i + 1) % 4, 'C-N') for i in range(4)])
    assert assembly_to_primary_sequence(a) == [['A', 'B', 'C', 'D']]


def test_repeated_identities_survive_but_unknowns_and_unassigned_are_excluded():
    a = assembly(['glycine', 'glycine', None, 'sugar'], [(0, 1, 'C-N'), (1, 2, 'C-C')])
    # The aggregate unassigned placeholder is not a sequencable MolNode.
    a.g.add_node(a.unassigned, molnode=None, tags={99}, identity=None)
    before = a.to_dict()
    readout = LinearReadout.from_assembly_graph(a)
    assert assembly_to_primary_sequence(a) == [['glycine', 'glycine'], ['sugar']]
    assert sum(map(len, readout.paths)) == 3
    assert 'X' not in readout.primary_sequence()
    assert a.to_dict() == before


def test_unidentified_bridge_splits_named_paths():
    a = assembly(['A', None, 'B'], [(0, 1, 'C-C'), (1, 2, 'C-N')])
    assert assembly_to_primary_sequence(a) == [['A'], ['B']]


def test_direct_and_legacy_readouts_split_at_unidentified_nodes():
    a = assembly(['A', None, 'B'], [(0, 1, 'C-C'), (1, 2, 'C-N')])
    nodes = a.monomer_nodes()
    direct = LinearReadout(a, [nodes])
    legacy = LinearReadout.from_dict({
        'assembly_graph': a.to_dict(), 'paths': [[node.to_dict() for node in nodes]],
    })
    for readout in (direct, legacy):
        assert readout.primary_sequence_parts() == [['A'], ['B']]
        assert all(node.is_identified for path in readout.paths for node in path)


def test_empty_assembly_and_unassigned_only_are_empty_sequences():
    a = assembly([])
    assert assembly_to_primary_sequence(a) == []
    a.g.add_node(a.unassigned, molnode=None, tags={1}, identity=None)
    assert assembly_to_primary_sequence(a) == []


def test_wholly_unidentified_assembly_is_an_empty_readout():
    a = assembly([None, None], [(0, 1, 'C-C')])
    readout = LinearReadout.from_assembly_graph(a)
    assert readout.paths == []
    assert readout.primary_sequence_parts() == []
    assert readout.primary_sequence() == []


@pytest.mark.parametrize('seed', range(10))
def test_path_partition_preserves_all_nodes_and_ignores_insertion_order(seed):
    topology = nx.gnp_random_graph(12, .25, seed=seed)
    a = assembly([f'unit-{i}' for i in topology], [(u, v, 'C-C') for u, v in topology.edges])
    readout = LinearReadout.from_assembly_graph(a)
    nodes = [node.enc for path in readout.paths for node in path]
    assert len(nodes) == len(set(nodes)) == 12
    by_enc = {data['molnode'].enc: n for n, data in a.g.nodes(data=True)}
    assert all(a.g.has_edge(by_enc[u.enc], by_enc[v.enc])
               for path in readout.paths for u, v in zip(path, path[1:]))

    rng = random.Random(seed)
    nodes = list(a.g.nodes(data=True))
    edges = list(a.g.edges(data=True))
    rng.shuffle(nodes)
    rng.shuffle(edges)
    shuffled = nx.Graph()
    shuffled.add_nodes_from(nodes)
    shuffled.add_edges_from((v, u, data) for u, v, data in edges)
    reordered = LinearReadout.from_assembly_graph(AssemblyGraph(shuffled))
    assert [[n.enc for n in path] for path in reordered.paths] == [[n.enc for n in path] for path in readout.paths]


@pytest.mark.parametrize('role', ['unknown', '', True, 1, []])
def test_invalid_sequence_role_is_rejected(role):
    with pytest.raises(ValueError, match='sequence_role'):
        MatchingRule('unit', 'C', {}, [], sequence_role=role)


def test_rule_role_defaults_roundtrips_and_changes_provenance(tmp_path):
    legacy = {'name': 'unit', 'smiles': 'C'}
    rule = MatchingRule.from_dict(legacy)
    assert rule.sequence_role is None
    tagged = replace(rule, sequence_role='tailoring')
    assert MatchingRule.from_dict(tagged.to_dict()).sequence_role == 'tailoring'
    assert RuleSet(False, [], [rule]).matching_rules_hash != RuleSet(False, [], [tagged]).matching_rules_hash
    matching = tmp_path / 'mxn.yml'
    matching.write_text('- name: unit\n  smiles: C\n  sequence_role: tailoring\n')
    reaction = tmp_path / 'rxn.yml'
    reaction.write_text('[]\n')
    loaded = RuleSet.load_from_files(reaction, matching)
    assert loaded.matching_rules[0].sequence_role == 'tailoring'


def test_default_modification_rules_have_tailoring_roles(ruleset):
    for name in ('glycosylation', 'methylation', 'oxidation', 'amination', 'phosphate', 'chlorination'):
        matching = [rule for rule in ruleset.matching_rules if rule.name == name]
        assert matching and all(rule.sequence_role == 'tailoring' for rule in matching)
    # A modified monomer remains a whole unit, not a separate tailoring event.
    assert all(rule.sequence_role is None for rule in ruleset.matching_rules if rule.name == 'N-methylalanine')


def test_sequences_from_sampled_assemblies_and_saved_result(ruleset):
    result = run_retromol(Submission('CCCCCCCCCC(=O)NCC(=O)NC(C)C(=O)O'), ruleset)
    assemblies = result.sample_assemblies(10, seed=42)
    sequences = [assembly_to_primary_sequence(a) for a in assemblies]
    assert sorted(len(parts[0]) for parts in sequences) == [3, 4, 5, 6, 7]
    assert all(len(parts) == 1 for parts in sequences)
    restored = Result.from_dict(json.loads(json.dumps(result.to_dict())))
    assert [assembly_to_primary_sequence(a) for a in restored.sample_assemblies(10, seed=42)] == sequences
