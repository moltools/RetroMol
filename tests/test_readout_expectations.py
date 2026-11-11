"""The broad integration fixtures must catch order, grouping and coverage errors."""

import pytest
import yaml

from .helpers_retromol import (
    ExpectedAssembly, assert_assemblies, assert_linear_readout,
    load_expected_readouts, parse_compound,
)
from .test_assembly_sampling import layered_result
from .test_readout import assembly


FINE = ExpectedAssembly(
    1.0, (("C", 1),) * 4,
    ((0, 1, ("C-C:SINGLE",)), (1, 2, ("C-C:SINGLE",)), (2, 3, ("C-C:SINGLE",))),
)
COARSE = ExpectedAssembly(1.0, (("CC", 2),) * 2, ((0, 1, ("C-C:SINGLE",)),))


def test_each_assembly_alternative_requires_its_own_readout():
    result, *_ = layered_result()
    expected = ((('C', 'C', 'C', 'C'),), (('CC', 'CC'),))
    assert_assemblies(result, [FINE, COARSE], expected_readouts=expected)
    with pytest.raises(AssertionError, match="Matching graphs with wrong readouts"):
        assert_assemblies(result, [FINE, COARSE], expected_readouts=(expected[0], (('CC',), ('CC',))))
    with pytest.raises(ValueError, match="one readout for every"):
        assert_assemblies(result, [FINE, COARSE], expected_readouts=expected[:1])


@pytest.mark.parametrize('wrong', [
    (('norleucine', 'alanine'),),
    (('alanine',), ('norleucine',)),
    (('alanine',),),
    (('alanine', 'norleucine', 'norleucine'),),
], ids=['reversed', 'split', 'missing', 'duplicated'])
def test_exact_readout_comparison_rejects_wrong_order_or_parts(wrong, ruleset):
    result = parse_compound('CC(N)C(=O)NC(CCCC)C(=O)O', ruleset)
    assert_linear_readout(result.assembly_graph, (('alanine', 'norleucine'),))
    with pytest.raises(AssertionError, match="Wrong linear readout"):
        assert_linear_readout(result.assembly_graph, wrong)


def test_identical_parts_are_not_deduplicated_and_unknown_units_are_omitted():
    graph = assembly(['repeat', 'repeat', None])
    assert_linear_readout(graph, (('repeat',), ('repeat',)))
    with pytest.raises(AssertionError, match="Wrong linear readout"):
        assert_linear_readout(graph, (('repeat',),))
    assert_linear_readout(assembly([None]), ())


@pytest.mark.parametrize('data, message', [
    ({}, 'exactly the assembly compound set'),
    ({'case': {'default_assembly': 0, 'assemblies': []}}, 'every expected assembly'),
    ({'case': {'default_assembly': 0, 'assemblies': [[['C', 'C', 'C']]]}}, 'every identified unit'),
    ({'case': {'default_assembly': 1, 'assemblies': [[['C', 'C', 'C', 'C']]]}}, 'invalid default assembly index'),
])
def test_readout_fixture_validation_rejects_incomplete_expectations(tmp_path, data, message):
    path = tmp_path / 'readouts.yml'
    path.write_text(yaml.safe_dump(data))
    with pytest.raises(ValueError, match=message):
        load_expected_readouts(path, {'case': (FINE,)})


def test_default_readout_can_describe_a_different_decomposition(tmp_path):
    path = tmp_path / 'readouts.yml'
    path.write_text(yaml.safe_dump({'case': {
        'default_sequence': [['C', 'C', 'C', 'C']],
        'assemblies': [[['CC', 'CC']]],
    }}))
    expected = load_expected_readouts(path, {'case': (COARSE,)})['case']
    assert expected.default == (('C', 'C', 'C', 'C'),)
    assert expected.assemblies == ((('CC', 'CC'),),)
