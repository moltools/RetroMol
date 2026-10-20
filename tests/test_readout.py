"""Linear assembly parts, role boundaries, and preservation of branched units."""

import json
import random
from dataclasses import replace

import networkx as nx
import pytest
from rdkit import Chem

from retromol.chem.mol import smiles_to_mol
from retromol.model.assembly_graph import AssemblyGraph, RootBondLink
from retromol.model.identity import MolIdentity
from retromol.model.reaction_graph import MolNode
from retromol.model.readout import (
    LinearReadout, _backbone_entry_tags, _connected_backbone_order,
    assembly_to_primary_sequence, split_named_paths,
)
from retromol.model.result import Result
from retromol.model.rules import MatchingRule, RuleSet
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol
from retromol_synthesis.reconstruction import reconstruct_linear_readout


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


@pytest.mark.parametrize('reverse_atoms', [False, True])
@pytest.mark.parametrize('smiles, expected', [
    ('NCC(=O)NC(C)C(=O)O', [['glycine', 'alanine']]),
    ('CC(N)C(=O)NCC(=O)O', [['alanine', 'glycine']]),
    ('NC(CC(=O)O)C(=O)NC(C)C(=O)O', [['aspartic acid', 'alanine']]),
    ('NCC(=O)NC(CC(=O)O)C(=O)O', [['glycine', 'aspartic acid']]),
    ('NCC(=O)N(C)C(C)C(=O)O', [['glycine', 'N-methylalanine']]),
    ('NCC(=O)NC(C)C(N)=O', [['glycine', 'alanine'], ['amination']]),
    ('CC(=O)CC(O)CC(=O)O', [['acetic acid starter', 'A1', 'B1']]),
    ('CCC(=O)CC(=O)O', [['propanoic acid starter', 'A1']]),
    ('NCC(=O)CC(=O)O', [['glycine', 'A1']]),
    ('CC(=O)CC(=O)NCC(=O)NC(C)C(=O)O', [['acetic acid starter', 'A1', 'glycine', 'alanine']]),
    ('CC(=O)CC(=O)NCC(=O)CC(=O)O', [['acetic acid starter', 'A1', 'glycine', 'A1']]),
])
def test_readout_runs_from_n_or_starter_end_to_carboxyl_end(smiles, expected, reverse_atoms, ruleset):
    mol = smiles_to_mol(smiles)
    if reverse_atoms:
        mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    result = run_retromol(Submission(Chem.MolToSmiles(mol, canonical=False)), ruleset)
    a = result.assembly_graph
    before = a.to_dict()
    readout = LinearReadout.from_assembly_graph(a)
    assert readout.primary_sequence_parts() == expected
    assert a.to_dict() == before
    assert result.calculate_coverage(a) == 1.0
    # Orientation must also work without the in-memory root molecule.
    restored = AssemblyGraph.from_dict(json.loads(json.dumps(before)))
    assert 'root_mol' not in restored.g.graph
    assert assembly_to_primary_sequence(restored) == expected
    # Direct construction and old saved readouts may contain the opposite direction.
    reversed_readout = readout.to_dict()
    reversed_readout['paths'] = [list(reversed(path)) for path in reversed_readout['paths']]
    assert LinearReadout.from_dict(reversed_readout).primary_sequence_parts() == expected
    assert LinearReadout(a, [list(reversed(path)) for path in readout.paths]).primary_sequence_parts() == expected


def test_oxygen_linked_backbone_uses_the_same_direction_rule(ruleset):
    backbone_rules = replace(ruleset, matching_rules=[
        replace(rule, sequence_role='backbone') if rule.name in ('glycine', 'lactic acid', 'alanine') else rule
        for rule in ruleset.matching_rules
    ])
    result = run_retromol(Submission('NCC(=O)OC(C)C(=O)NC(C)C(=O)O'), backbone_rules)
    assert result.linear_readout.primary_sequence_parts() == [['glycine', 'lactic acid', 'alanine']]


@pytest.mark.parametrize('smiles, expected', [
    ('NCC(=O)NC(COC(=O)C(C)N)C(=O)O', [['glycine', 'serine'], ['alanine']]),
    ('NCC(=O)N(C)C(COC(=O)C(C)N)C(=O)O', [['glycine', 'N-methylserine'], ['alanine']]),
    ('NC(COC(=O)C(C)N)C(=O)O', [['alanine', 'serine']]),
])
def test_amino_entry_connection_distinguishes_side_esters_from_continuations(smiles, expected, ruleset):
    # All units are explicitly backbone, so separation cannot come from the
    # fallback bond filter. An unconnected amine alone does not exclude an ester.
    backbone_rules = replace(ruleset, matching_rules=[
        replace(rule, sequence_role='backbone')
        if rule.name in ('glycine', 'alanine', 'serine', 'N-methylserine') else rule
        for rule in ruleset.matching_rules
    ])
    result = run_retromol(Submission(smiles), backbone_rules)
    assert result.linear_readout.primary_sequence_parts() == expected
    assert result.calculate_coverage() == 1.0


def test_amino_entry_excludes_lysine_side_chain_nitrogen():
    tagged = '[1NH2][2CH]([3CH2][4CH2][5CH2][6CH2][7NH2])[8C](=[9O])[10OH]'
    node = MolNode('lysine', smiles_to_mol(tagged), 'NC(CCCCN)C(=O)O')
    assert _backbone_entry_tags(node) == {1}


def test_retained_terminal_carboxyl_anchors_alternative_pk_decompositions(ruleset):
    # This structure also admits a reverse-cut route containing malonic acid.
    # Its competing artificial caps must not put the retained terminal acid first.
    result = run_retromol(Submission('CCC(=O)CC(O)CC(=O)O'), ruleset)
    sequences = [assembly_to_primary_sequence(a) for a in result.sample_assemblies(10, seed=42)]
    assert sorted(sequences) == sorted([
        [['propanoic acid starter', 'A1', 'B1']],
        [['propanoic acid starter', 'B1', 'malonic acid']],
    ])


def test_carboxyl_direction_does_not_reorient_tailoring_paths(ruleset):
    tailoring_rules = replace(ruleset, matching_rules=[
        replace(rule, sequence_role='tailoring') for rule in ruleset.matching_rules
    ])
    result = run_retromol(Submission('NCC(=O)NC(C)C(=O)O'), tailoring_rules)
    assert result.linear_readout.primary_sequence_parts() == [['alanine', 'glycine']]


@pytest.mark.parametrize('explicit_roles', [False, True])
@pytest.mark.parametrize('smiles, expected', [
    # Gly enters the backbone N; Ala acylates the side-chain N.
    ('NCC(=O)NC(CCCCNC(=O)C(N)C)C(=O)O', [['glycine', 'lysine'], ['alanine']]),
    # Exchanging the acyl donors must exchange their roles in the readout.
    ('CC(N)C(=O)NC(CCCCNC(=O)CN)C(=O)O', [['alanine', 'lysine'], ['glycine']]),
    # A side-chain amide alone is insufficient evidence to split a path.
    ('NC(CCCCNC(=O)C(N)C)C(=O)O', [['alanine', 'lysine']]),
])
def test_side_chain_amides_compete_with_connected_backbone_entries(smiles, expected, explicit_roles, ruleset):
    if explicit_roles:
        ruleset = replace(ruleset, matching_rules=[
            replace(rule, sequence_role='backbone') if rule.name in ('glycine', 'alanine', 'lysine') else rule
            for rule in ruleset.matching_rules
        ])
    result = run_retromol(Submission(smiles), ruleset)
    before = result.assembly_graph.to_dict()
    assert result.linear_readout.primary_sequence_parts() == expected
    assert result.assembly_graph.to_dict() == before
    restored = AssemblyGraph.from_dict(json.loads(json.dumps(before)))
    assert assembly_to_primary_sequence(restored) == expected


def test_cyclic_peptide_preserves_conventional_start_and_carboxyl_direction(ruleset):
    # Around this ring: Gly -> Ala -> Ser -> Gly; there is no natural terminus.
    result = run_retromol(Submission('O=C1NCC(=O)NC(C)C(=O)NC(CO)1'), ruleset)
    readout = result.linear_readout
    assert readout.primary_sequence_parts() == [['alanine', 'serine', 'glycine']]
    restored = LinearReadout.from_dict(readout.to_dict())
    assert restored.primary_sequence_parts() == readout.primary_sequence_parts()


def crosslinked_capped_chain():
    """A source, two capped units and an amino acid, with competing ring bonds."""
    g = nx.Graph()
    for index, (name, smiles) in enumerate([
        ('z-source', '[1CH3]C(=O)O'),
        ('b-unit', '[2CH2](SO)[3CH2]C(=O)O'),
        ('a-unit', '[4CH2](SO)[5CH2]C(=O)O'),
        ('c-end', '[6NH2][7CH]([8CH3])[9C](=[10O])[11OH]'),
    ]):
        mol = smiles_to_mol(smiles)
        identity = MolIdentity(MatchingRule(name, smiles, {}, [], sequence_role='backbone'))
        node = MolNode(f'unit-{index}', mol, smiles, identity, True)
        g.add_node(index, molnode=node, identity=identity,
                   tags={atom.GetIsotope() for atom in mol.GetAtoms() if atom.GetIsotope()})
    for u, v, a, b, symbol in [
        (0, 1, 1, 2, 'C'), (1, 2, 3, 4, 'C'), (2, 3, 5, 6, 'N'),
        (3, 1, 9, 2, 'C'),  # Ring closure into an already connected entry.
        (1, 3, 3, 6, 'N'),  # Exit-to-entry shortcut that skips a unit.
    ]:
        bond = RootBondLink(a - 1, b - 1, a, b, 'C', symbol, 'SINGLE', 1.0)
        if g.has_edge(u, v):
            g[u][v]['bonds'].append(bond)
            g[u][v]['n_bonds'] += 1
        else:
            g.add_edge(u, v, bonds=[bond], n_bonds=1)
    return AssemblyGraph(g)


def test_entry_exit_connections_select_order_before_undirected_crosslinks():
    a = crosslinked_capped_chain()
    before = a.to_dict()
    expected = [['z-source', 'b-unit', 'a-unit', 'c-end']]
    assert assembly_to_primary_sequence(a) == expected
    # This connected path takes a ring shortcut; neither orientation is correct.
    saved = LinearReadout(a, [[a.g.nodes[i]['molnode'] for i in (0, 1, 3, 2)]])
    assert saved.primary_sequence_parts() == expected
    assert LinearReadout.from_dict(saved.to_dict()).primary_sequence_parts() == expected
    assert a.to_dict() == before

    # Labels and graph insertion order cannot choose a different backbone.
    for i, data in a.g.nodes(data=True):
        rule = replace(data['molnode'].identity.matched_rule, name=f'renamed-{3 - i}')
        data['identity'] = MolIdentity(rule)
        data['molnode'] = replace(data['molnode'], identity=data['identity'])
    reordered = nx.Graph()
    reordered.add_nodes_from(reversed(list(a.g.nodes(data=True))))
    reordered.add_edges_from((v, u, data) for u, v, data in reversed(list(a.g.edges(data=True))))
    readout = LinearReadout.from_assembly_graph(AssemblyGraph(reordered))
    assert [[node.enc for node in path] for path in readout.paths] == [
        ['unit-0', 'unit-1', 'unit-2', 'unit-3'],
    ]


def test_competing_complete_port_orders_are_left_ambiguous():
    a = crosslinked_capped_chain()
    # The terminal acid can also enter the middle unit: both 0-1-2-3 and
    # 0-1-3-2 now follow exit-to-entry links, so neither is uniquely supported.
    a.g[3][2]['bonds'].append(RootBondLink(8, 3, 9, 4, 'C', 'C', 'SINGLE', 1.0))
    a.g[3][2]['n_bonds'] += 1
    assert _connected_backbone_order(a.g) is None


def test_reconstruction_consumes_the_oriented_peptide_readout(ruleset):
    smiles = 'NCC(=O)NC(C)C(=O)O'
    result = run_retromol(Submission(smiles), ruleset)
    reconstruction, = reconstruct_linear_readout(result)
    assert [name for name, _ in reconstruction.primary_sequence] == ['glycine', 'alanine']
    mol = smiles_to_mol(reconstruction.tagged_backbone_smiles)
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    assert Chem.MolToSmiles(mol) == Chem.MolToSmiles(smiles_to_mol(smiles))


ARTIFICIAL_ERYTHROMYCIN = (
    'CC[C@@H]1[C@@](C)(O)[C@H](O)[C@@H](C)C([C@H](C)C[C@@](C)(OC(C)=O)'
    '[C@H](O[C@H]2[C@H](O)[C@@H](N(C)C)C[C@@H](C)O2)[C@@H](C)'
    '[C@H](O[C@H]3C[C@@](C)(O)[C@@H](O)[C@H](C)O3)[C@@H](C)C(O1)=O)=O'
)


@pytest.mark.parametrize('stereo', [False, True])
@pytest.mark.parametrize('reverse_atoms', [False, True])
def test_macrolactone_closure_does_not_shortcut_the_linear_backbone(stereo, reverse_atoms, ruleset, stereo_ruleset):
    mol = smiles_to_mol(ARTIFICIAL_ERYTHROMYCIN)
    if reverse_atoms:
        mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    result = run_retromol(Submission(Chem.MolToSmiles(mol, canonical=False)), stereo_ruleset if stereo else ruleset)
    a = result.assembly_graph
    before = a.to_dict()
    backbone = ['propanoic acid starter', 'B^R6^S', 'B^S2^S', 'A2^R', 'D6^S', 'B^S2^S', 'B^S2^S'] if stereo else [
        'propanoic acid starter', 'B6', 'B2', 'A2', 'D6', 'B2', 'B2',
    ]
    expected = [backbone, ['acetic acid'], ['glycosylation'], ['glycosylation']]
    readout = LinearReadout.from_assembly_graph(a)
    assert readout.primary_sequence_parts() == expected
    assert result.calculate_coverage(a) == 1.0
    assert a.to_dict() == before
    assert assembly_to_primary_sequence(AssemblyGraph.from_dict(json.loads(json.dumps(before)))) == expected
    # Every sequential connection follows C-C; the original lactone closure is
    # still present as a C-O connection between the first and last PK units.
    path = readout.paths[0]
    assert a.g.has_edge(path[1].enc, path[-1].enc)
    for left, right in zip(path, path[1:]):
        assert all((b.a1_symbol, b.a2_symbol) == ('C', 'C') for b in a.g[left.enc][right.enc]['bonds'])

    # Direction/topology must not depend on a name containing "starter".
    renamed = a.g.copy()
    for i, node in enumerate(path):
        identity = MolIdentity(replace(node.identity.matched_rule, name=f'unit-{len(path) - i}'))
        renamed.nodes[node.enc]['molnode'] = replace(node, identity=identity)
        renamed.nodes[node.enc]['identity'] = identity
    renamed_path = LinearReadout.from_assembly_graph(AssemblyGraph(renamed)).paths[0]
    assert [node.enc for node in renamed_path] == [node.enc for node in path]

    # A side attachment must not become a continuation merely because it is
    # the only link between two chains. The acceptor still has its own entry.
    opened = a.g.copy()
    opened.remove_edge(path[3].enc, path[4].enc)
    opened_readout = LinearReadout.from_assembly_graph(AssemblyGraph(opened))
    closure_pair = {path[1].enc, path[-1].enc}
    assert not any({left.enc, right.enc} == closure_pair
                   for part in opened_readout.paths for left, right in zip(part, part[1:]))
    assert opened_readout.primary_sequence_parts() == [backbone[:4], backbone[4:], *expected[1:]]
