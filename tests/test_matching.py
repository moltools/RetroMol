"""Whole-fragment matching with optional conditions on the tagged structure."""

import hashlib
import json
from dataclasses import replace

import pytest
import yaml

from retromol.chem.matching import match_mol
from retromol.chem.mol import smiles_to_mol
from retromol.model.result import Result
from retromol.model.rules import MatchingCase, MatchingRule, RuleSet
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol
from retromol_antismash.modules import NRPSModule, module_primary_sequence_tokens
from retromol_synthesis.reconstruction import reconstruct_named_sequence


ADDED_ACID = '[0C](=[0O])[0OH]'


def starter_rule():
    return MatchingRule('acetic acid starter', 'CC(=O)O', {}, [],
                        requires_smarts=ADDED_ACID, sequence_role='backbone')


def acid_rule():
    return MatchingRule.from_dict({
        'name': 'acetic acid', 'smiles': 'CC(=O)O', 'sequence_role': 'tailoring',
        'cases': [{'requires_smarts': ADDED_ACID, 'name': 'acetic acid starter',
                   'sequence_role': 'backbone'}],
    })


def test_single_structural_rule_resolves_each_fragment_without_mutating_defaults():
    rule = acid_rule()
    original_definition = rule.to_dict()
    unrelated = MatchingRule('methane', 'C', {}, [])
    for rules in ([rule, unrelated], [unrelated, rule]):
        for fragment, name, role in [
            (smiles_to_mol('[1CH3]C(=O)O'), 'acetic acid starter', 'backbone'),
            (smiles_to_mol('[1CH3][2C](=[3O])O'), 'acetic acid', 'tailoring'),
            (Submission('CC(=O)O').mol, 'acetic acid', 'tailoring'),
        ]:
            assert rule.is_match(fragment)
            identity = match_mol(fragment, rules)
            assert identity.name == name
            assert identity.matched_rule.sequence_role == role
            assert identity.matched_rule.cases == []
    assert rule.to_dict() == original_definition
    # Cases are evaluated only after the whole structure and any top-level gate match.
    assert rule.match(smiles_to_mol('[1CH3][2CH2]C(=O)O')) is None
    assert replace(rule, requires_smarts='N').match(smiles_to_mol('[1CH3]C(=O)O')) is None


def test_cases_use_first_matching_override_and_inherit_unspecified_defaults():
    rule = replace(acid_rule(), cases=[
        MatchingCase('[0C]', name='first'),
        MatchingCase('[0O]', name='second', sequence_role='backbone'),
    ])
    generated = smiles_to_mol('[1CH3]C(=O)O')
    resolved = rule.match(generated)
    assert (resolved.name, resolved.sequence_role) == ('first', 'tailoring')
    role_only = replace(rule, cases=[MatchingCase(ADDED_ACID, sequence_role='backbone')])
    resolved = role_only.match(generated)
    assert (resolved.name, resolved.sequence_role) == ('acetic acid', 'backbone')


def test_case_queries_use_the_requested_chirality():
    rule = MatchingRule('alanine', 'NC(C)C(=O)O', {}, [], cases=[
        MatchingCase('N[C@@H](C)C(=O)O', name='stereo case'),
    ])
    same = smiles_to_mol('N[C@@H](C)C(=O)O')
    opposite = smiles_to_mol('N[C@H](C)C(=O)O')
    assert rule.match(same, match_stereochemistry=True).name == 'stereo case'
    assert rule.match(opposite, match_stereochemistry=True).name == 'alanine'
    assert rule.match(opposite, match_stereochemistry=False).name == 'stereo case'


@pytest.mark.parametrize('cases', [
    None, {}, 'case', [None], [{}], [{'requires_smarts': ''}],
    [{'requires_smarts': '[C'}], [{'requires_smarts': True}],
    [{'requires_smarts': 'C', 'smiles': 'C'}],
    [{'requires_smarts': 'C', 'name': ''}],
    [{'requires_smarts': 'C', 'sequence_role': 'unknown'}],
])
def test_invalid_cases_fail_when_loading(cases):
    with pytest.raises(ValueError, match='case'):
        MatchingRule.from_dict({'name': 'invalid', 'smiles': 'C', 'cases': cases})


def test_cases_roundtrip_and_catalog_names(tmp_path):
    rule = acid_rule()
    restored = MatchingRule.from_dict(json.loads(json.dumps(rule.to_dict())))
    assert restored.to_dict() == rule.to_dict()
    assert restored.id == replace(rule, cases=[]).id
    ruleset = RuleSet(False, [], [restored])
    assert ruleset.matching_rules_hash != RuleSet(False, [], [replace(rule, cases=[])]).matching_rules_hash
    assert [(r.name, r.sequence_role) for r in ruleset.matching_rule_variants] == [
        ('acetic acid', 'tailoring'), ('acetic acid starter', 'backbone'),
    ]
    assert [r.name for r in ruleset.find_structural_matches(smiles_to_mol('[1CH3]C(=O)O'))] == [
        'acetic acid starter',
    ]
    assert [r.name for r in ruleset.find_structural_matches(Submission('CC(=O)O').mol)] == [
        'acetic acid',
    ]
    reaction_path = tmp_path / 'rxn.yml'
    matching_path = tmp_path / 'mxn.yml'
    reaction_path.write_text('[]\n')
    matching_path.write_text(yaml.safe_dump([rule.to_dict()]))
    loaded = RuleSet.load_from_files(reaction_path, matching_path)
    assert len(loaded.matching_rules) == 1
    assert loaded.matching_rules_hash == ruleset.matching_rules_hash
    assert match_mol(smiles_to_mol('[1CH3]C(=O)O'), loaded.matching_rules).name == 'acetic acid starter'


@pytest.mark.parametrize('name', ['acetic acid', 'propanoic acid'])
def test_default_starter_cases_remain_available_for_named_reconstruction(name, ruleset):
    assert [r.name for r in ruleset.matching_rules if r.name in (name, name + ' starter')] == [name]
    reconstructed = reconstruct_named_sequence(ruleset, [name + ' starter', 'A1'])
    assert reconstructed.tagged_backbone_smiles
    assert [part[0] for part in reconstructed.primary_sequence] == [name + ' starter', 'A1']


@pytest.mark.parametrize('name, prefix', [
    ('acetic acid', '[1CH3]'), ('propanoic acid', '[1CH3][2CH2]'),
])
@pytest.mark.parametrize('acid, starter', [
    ('C(=O)O', True),  # Complete cap added by PK cleavage.
    ('[3C](=O)O', True),  # Original thiazole carbon; both oxygens restored.
    ('[3C](=[4O])O', False),  # O-/N-acyl cleavage retains carbonyl oxygen.
    ('[3C](=O)[5OH]', False),  # One new oxygen alone is insufficient.
    ('[3C](=[4O])[5OH]', False),  # Original free acid.
])
def test_default_acid_case_accepts_restored_oxygens_with_either_carbon_origin(name, prefix, acid, starter, ruleset):
    rule, = [r for r in ruleset.matching_rules if r.name == name]
    identity = match_mol(smiles_to_mol(prefix + acid), [rule])
    assert identity.name == name + (' starter' if starter else '')
    assert identity.matched_rule.sequence_role == ('backbone' if starter else 'tailoring')


@pytest.mark.parametrize('smiles, expected', [
    ('[1CH3]C(=O)O', True),
    ('[1CH3][2C](=[3O])O', False),  # Cleaved acetyl group; only OH is added.
    ('[1CH3][2C](=[3O])[4OH]', False),  # Original free acid.
    ('[1CH3]C(=[3O])O', False),  # A partially artificial cap is insufficient.
    ('[1CH3]C(=O)[4OH]', False),
    ('[1CH3][2CH2]C(=O)O', False),  # Query matches, whole-fragment SMILES does not.
])
def test_requires_smarts_checks_tagged_atoms_without_replacing_structure(smiles, expected):
    assert starter_rule().is_match(smiles_to_mol(smiles)) is expected


def test_conditional_rule_precedes_ordinary_acid_fallback():
    starter = starter_rule()
    ordinary = MatchingRule('acetic acid', 'CC(=O)O', {}, [], sequence_role='tailoring')
    generated = smiles_to_mol('[1CH3]C(=O)O')
    original = Submission('CC(=O)O').mol
    assert match_mol(generated, [starter, ordinary]).name == 'acetic acid starter'
    assert match_mol(original, [starter, ordinary]).name == 'acetic acid'
    # The feature does not change first-match ordering.
    assert match_mol(generated, [ordinary, starter]).name == 'acetic acid'


def test_requires_smarts_uses_the_requested_chirality():
    rule = MatchingRule('conditional alanine', 'NC(C)C(=O)O', {}, [],
                        requires_smarts='N[C@@H](C)C(=O)O')
    same = smiles_to_mol('N[C@@H](C)C(=O)O')
    opposite = smiles_to_mol('N[C@H](C)C(=O)O')
    assert rule.is_match(same, match_stereochemistry=True)
    assert not rule.is_match(opposite, match_stereochemistry=True)
    assert rule.is_match(opposite, match_stereochemistry=False)


@pytest.mark.parametrize('condition', ['', '   ', True, 0, [], {}, '[C'])
def test_invalid_requires_smarts_fails_when_loading(condition):
    with pytest.raises(ValueError, match='requires_smarts'):
        MatchingRule.from_dict({'name': 'invalid', 'smiles': 'C', 'requires_smarts': condition})


def test_condition_serialization_ids_and_legacy_defaults(tmp_path):
    starter = starter_rule()
    restored = MatchingRule.from_dict(json.loads(json.dumps(starter.to_dict())))
    assert restored.requires_smarts == ADDED_ACID
    assert restored.id == starter.id
    assert restored.is_match(smiles_to_mol('[1CH3]C(=O)O'))
    assert not restored.is_match(Submission('CC(=O)O').mol)

    legacy = MatchingRule.from_dict({'name': 'acetic acid', 'smiles': 'CC(=O)O'})
    assert legacy.requires_smarts is None
    assert legacy.is_match(Submission('CC(=O)O').mol)
    assert legacy.id == hashlib.sha256(b'CC(=O)O').hexdigest()
    assert legacy.id != starter.id
    assert replace(starter, requires_smarts='[0C]').id != starter.id
    assert RuleSet(False, [], [starter]).matching_rules_hash != RuleSet(
        False, [], [replace(starter, requires_smarts=None)],
    ).matching_rules_hash

    reaction_path = tmp_path / 'rxn.yml'
    matching_path = tmp_path / 'mxn.yml'
    reaction_path.write_text('[]\n')
    matching_path.write_text(yaml.safe_dump([starter.to_dict(), legacy.to_dict()]))
    loaded = RuleSet.load_from_files(reaction_path, matching_path)
    assert match_mol(smiles_to_mol('[1CH3]C(=O)O'), loaded.matching_rules).name == starter.name
    assert match_mol(Submission('CC(=O)O').mol, loaded.matching_rules).name == legacy.name


@pytest.mark.parametrize('smiles, expected', [
    ('CC(=O)CC(=O)O', [['acetic acid starter', 'A1']]),
    ('CCC(=O)CC(=O)O', [['propanoic acid starter', 'A1']]),
    ('CC(O)CC(=O)O', [['acetic acid starter', 'B1']]),
    ('CC=CC(=O)O', [['acetic acid starter', 'C1']]),
    ('CCCC(=O)O', [['acetic acid starter', 'D1']]),
    ('CC(=O)CC(=O)CC(=O)O', [['acetic acid starter', 'A1', 'A1']]),
    ('CC(OC(C)=O)CC(=O)O', [['acetic acid starter', 'B1'], ['acetic acid']]),
    ('NC(COC(C)=O)C(=O)NC(C)C(=O)O', [['serine', 'alanine'], ['acetic acid']]),
    ('NC(CCCCNC(C)=O)C(=O)NC(C)C(=O)O', [['lysine', 'alanine'], ['acetic acid']]),
    ('CC(=O)NCC(=O)NC(C)C(=O)O', [['glycine', 'alanine'], ['acetic acid']]),
    ('CC(=O)O', [['acetic acid']]),
    ('CCC(=O)O', [['propanoic acid']]),
])
def test_default_rules_distinguish_pk_starters_and_acyl_tailoring(smiles, expected, ruleset):
    result = run_retromol(Submission(smiles), ruleset)
    assert result.linear_readout.primary_sequence_parts() == expected
    assert result.calculate_coverage() == 1.0
    restored = Result.from_dict(json.loads(json.dumps(result.to_dict())))
    assert restored.linear_readout.primary_sequence_parts() == expected
    assert restored.calculate_coverage() == 1.0


def test_old_results_keep_their_saved_identities(ruleset):
    # Before conditional starter rules, this same PK fragment was ordinary acid.
    old_rules = replace(ruleset, matching_rules=[replace(r, cases=[]) for r in ruleset.matching_rules])
    result = run_retromol(Submission('CC(=O)CC(=O)O'), old_rules)
    data = result.to_dict()
    for identity in data['reaction_graph']['identities']:
        identity.pop('requires_smarts')
        identity.pop('cases')
    restored = Result.from_dict(data)
    assert restored.linear_readout.primary_sequence_parts() == [['A1'], ['acetic acid']]


@pytest.mark.parametrize('smiles, name', [
    ('CC(=O)O', 'acetic acid'),
    ('CCC(=O)O', 'propanoic acid'),
])
def test_predicted_substrate_smiles_do_not_imply_generated_atoms(smiles, name, ruleset):
    module = NRPSModule.from_dict({
        'module_index_in_gene': 0, 'start': 0, 'end': 100,
        'gene_id': 'test', 'gene_strand': '+', 'present_domains': [],
        'predicted_substrate': {'smiles': smiles},
    })
    assert module_primary_sequence_tokens(module, ruleset) == (name, [name])
