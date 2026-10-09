"""Single-mode drawing limits use distinct, reproducible assembly projections."""

import json
from pathlib import Path

import pytest

from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.assembly_sampling import AssemblyLimitError, AssemblySpace
from retromol.model.result import Result
from retromol.visualization.assembly_layout import MoleculeLayout
from retromol_cli import cli as cli_module


SMILES = 'CCCCCCCCCC(=O)NCC(=O)NC(C)C(=O)O'


@pytest.fixture
def single_run(monkeypatch, tmp_path):
    # Exercise parsing, serialization, selection and coverage; only the expensive
    # rendering and global logging configuration are replaced.
    drawn = []
    visualized = []
    layouts = []
    monkeypatch.setattr(cli_module, 'setup_logging', lambda **kwargs: None)
    monkeypatch.setattr(cli_module, 'add_file_handler', lambda *args, **kwargs: None)

    def record_drawing(assembly, **kwargs):
        drawn.append((Path(kwargs['savepath']).name, tuple(sorted(assembly.monomer_ids()))))
        assert kwargs['show_unassigned']
        layout = kwargs['layout']
        assert isinstance(layout, MoleculeLayout)
        if layouts:
            assert layout is layouts[0]
        layouts.append(layout)

    monkeypatch.setattr(AssemblyGraph, 'draw', record_drawing)
    monkeypatch.setattr(cli_module, 'visualize_reaction_graph', lambda *args, **kwargs: visualized.append(kwargs))

    def run(*options, smiles=SMILES):
        drawn.clear()
        visualized.clear()
        layouts.clear()
        output = tmp_path / 'output'
        monkeypatch.setattr('sys.argv', ['retromol', '-o', str(output), 'single', '-s', smiles, *options])
        cli_module.main()
        payload = json.loads((output / 'result.json').read_text())
        assert 'linear_readout' not in payload
        assert 'assembly_graph' not in payload
        assert len(visualized) == 1
        return Result.from_dict(payload), list(drawn)

    return run


def test_single_default_keeps_default_assembly_filename_and_selection(single_run):
    result, drawn = single_run()
    assert drawn == [('assembly_graph.png', tuple(sorted(result.assembly_graph.monomer_ids())))]


@pytest.mark.parametrize('maximum,expected_count', [(0, 0), (1, 1), (3, 3), (10, 5)])
def test_single_respects_drawing_limit_and_available_assemblies(single_run, maximum, expected_count):
    result, drawn = single_run('--max-assembly-graphs', str(maximum))
    assert len(drawn) == expected_count
    assert len({nodes for _, nodes in drawn}) == expected_count
    if maximum > 1:
        assert [name for name, _ in drawn] == [f'assembly_graph_{index:03d}.png' for index in range(1, expected_count + 1)]
        expected = [tuple(sorted(a.monomer_ids())) for a in result.sample_assemblies(maximum, seed=42)]
        assert [nodes for _, nodes in drawn] == expected


def test_single_passes_seed_to_sampling(single_run):
    result, first = single_run('--max-assembly-graphs', '3', '--seed', '17')
    _, repeated = single_run('--max-assembly-graphs', '3', '--seed', '17')
    assert first == repeated
    assert [nodes for _, nodes in first] == [tuple(sorted(a.monomer_ids())) for a in result.sample_assemblies(3, seed=17)]


def test_single_rerun_removes_stale_generated_drawings(single_run, tmp_path):
    output = tmp_path / 'output'
    output.mkdir()
    stale = [output / 'assembly_graph.png', *(output / f'assembly_graph_{i:03d}.png' for i in range(6, 11))]
    for path in stale:
        path.write_bytes(b'old drawing')
    user_image = output / 'assembly_graph_notes.png'
    user_image.write_bytes(b'user image')
    _, drawn = single_run('--max-assembly-graphs', '10')
    assert len(drawn) == 5
    assert not any(path.exists() for path in stale)
    assert user_image.read_bytes() == b'user image'


@pytest.mark.parametrize('maximum', ['-1', '1.5', 'invalid'])
def test_single_rejects_invalid_limit_before_creating_output(monkeypatch, tmp_path, capsys, maximum):
    output = tmp_path / 'output'
    monkeypatch.setattr('sys.argv', ['retromol', '-o', str(output), 'single', '-s', SMILES, '--max-assembly-graphs', maximum])
    with pytest.raises(SystemExit) as error:
        cli_module.main()
    assert error.value.code == 2
    assert 'must be a nonnegative integer' in capsys.readouterr().err
    assert not output.exists()


def test_single_enumeration_limit_fails_clearly_and_preserves_result(single_run, monkeypatch, tmp_path, caplog):
    def over_budget(*args, **kwargs):
        raise AssemblyLimitError('enumeration budget exceeded')

    monkeypatch.setattr(AssemblySpace, 'sample', over_budget)
    with pytest.raises(SystemExit) as error:
        single_run('--max-assembly-graphs', '10')
    assert error.value.code == 1
    assert 'enumeration budget exceeded' in caplog.text
    assert 'reaction graph is saved' in caplog.text
    assert Result.from_dict(json.loads((tmp_path / 'output' / 'result.json').read_text())).root


def test_daptomycin_samples_only_full_coverage_and_keeps_all_resolutions(single_run, caplog):
    from .data.integration_retromol import CASES

    smiles = next(smiles for name, smiles, *_ in CASES if name == 'daptomycin')
    with caplog.at_level('INFO'):
        result, drawn = single_run('--max-assembly-graphs', '20', smiles=smiles)
    space = result.assembly_space()
    assert space.count() == 5
    assert space.count(min_coverage=None) == len(drawn) == 5
    assert 'Drawing 5 of 5 eligible assemblies (5 total' in caplog.text
    assemblies = [space.build(frozenset(nodes)) for _, nodes in drawn]
    assert all(result.calculate_coverage(a) == 1 for a in assemblies)
    fatty_acids = {'decanoic acid', 'octanoic acid', 'hexanoic acid', 'butanoic acid', 'acetic acid'}
    assert fatty_acids <= {node.identity.name for a in assemblies for node in a.monomer_nodes()}
    assert len({nodes for _, nodes in drawn}) == 5
    threonine = next(n.enc for n in result.reaction_graph.identified_nodes.values()
                     if n.identity.name == 'threonine')
    for assembly in space.sample(20, seed=42, min_coverage=0):
        assert threonine in assembly.monomer_ids()
        assert 'B12' not in {n.identity.name for n in assembly.monomer_nodes()}
        assert result.calculate_coverage(assembly) == 1
    # The source graph remains agnostic, including the bypass's identified B12.
    assert any(n.identity.name == 'B12' for n in result.reaction_graph.identified_nodes.values())


def test_single_default_also_uses_maximum_if_default_projection_is_lower(single_run, monkeypatch):
    from .test_assembly_sampling import partially_identified_result

    parsed = partially_identified_result()
    monkeypatch.setattr(cli_module, 'run_retromol_with_timeout', lambda *args: parsed)
    result, drawn = single_run()
    assert result.calculate_coverage() == .75
    assert len(drawn) == 1
    assert result.calculate_coverage(result.assembly_space().build(frozenset(drawn[0][1]))) == 1
    _, lower = single_run('--min-coverage', '.75')
    assert lower == [('assembly_graph.png', tuple(sorted(result.assembly_graph.monomer_ids())))]
    _, all_views = single_run('--max-assembly-graphs', '10', '--min-coverage', '.75')
    assert len(all_views) == 4


@pytest.mark.parametrize('maximum', [1, 10])
def test_single_unattainable_minimum_draws_nothing_and_keeps_result(single_run, caplog, maximum):
    result, drawn = single_run('--max-assembly-graphs', str(maximum), '--min-coverage', '1',
                               smiles=SMILES + '.[Na+]')
    assert not drawn
    assert result.assembly_space().max_coverage < 1
    assert 'No assemblies meet minimum coverage 100.00%' in caplog.text


@pytest.mark.parametrize('minimum', ['-0.1', '1.1', 'NaN', 'inf', 'invalid'])
def test_single_rejects_invalid_coverage_before_creating_output(monkeypatch, tmp_path, capsys, minimum):
    output = tmp_path / 'output'
    monkeypatch.setattr('sys.argv', ['retromol', '-o', str(output), 'single', '-s', SMILES, '--min-coverage', minimum])
    with pytest.raises(SystemExit) as error:
        cli_module.main()
    assert error.value.code == 2
    assert 'must be a finite fraction between 0 and 1' in capsys.readouterr().err
    assert not output.exists()
