"""Molecular centroids keep alternative assembly drawings spatially comparable."""

import json

import matplotlib.pyplot as plt
import pytest
from rdkit import Chem

from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.result import Result
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol
from retromol.visualization.assembly_layout import MoleculeLayout

from .test_assembly_sampling import layered_result
from .test_result import make_result


def test_positions_average_original_heavy_atoms_and_place_real_unassigned_regions():
    result = make_result('[2H]CCO', [('[1H][2CH2][3CH2][99CH2]O', True)])
    assembly = result.assembly_graph
    layout = MoleculeLayout.from_mol(result.root.mol)
    assert set(layout.atom_positions) == {2, 3, 4}
    positions = assembly.monomer_positions(layout)
    expected = tuple((layout.atom_positions[2][axis] + layout.atom_positions[3][axis]) / 2 for axis in [0, 1])
    assert positions[assembly.monomer_ids()[0]] == pytest.approx(expected)
    assert positions[assembly.unassigned] == layout.atom_positions[4]
    assert layout.centroid([1, 99, 0]) is None
    assert layout.centroid([2, 2, 3, 99]) == pytest.approx(expected)


def test_alternatives_keep_common_nodes_and_plot_bounds_at_identical_coordinates():
    result, left, right, singles = layered_result()
    space = result.assembly_space()
    coarse = space.build(frozenset([left, right]))
    fine = space.build(frozenset([left, *singles[2:]]))
    layout = MoleculeLayout.from_mol(result.root.mol)
    assert coarse.monomer_positions(layout)[left] == fine.monomer_positions(layout)[left]
    assert coarse.unassigned not in coarse.monomer_positions(layout)
    fig, axes = plt.subplots(1, 2)
    try:
        coarse.draw(layout=layout, ax=axes[0], show_unassigned=True)
        fine.draw(layout=layout, ax=axes[1], show_unassigned=True)
        for ax in axes:
            assert ax.get_xlim() == layout.bounds[:2]
            assert ax.get_ylim() == layout.bounds[2:]
            assert ax.get_aspect() == 1
        assert len(axes[0].collections[0].get_offsets()) == 2
        assert len(axes[1].collections[0].get_offsets()) == 3
    finally:
        plt.close(fig)


def test_depiction_is_lazy_and_never_changes_source_or_saved_payload(monkeypatch):
    result, *_ = layered_result()
    before = json.loads(json.dumps(result.to_dict()))
    conformers = result.root.mol.GetNumConformers()
    with monkeypatch.context() as context:
        context.setattr(MoleculeLayout, 'from_mol', lambda *args: pytest.fail('Unexpected depiction during parsing/readout'))
        assemblies = result.sample_assemblies(10)
        _ = result.linear_readout
        assert json.loads(json.dumps(result.to_dict())) == before
    layout = assemblies[0].molecular_layout()
    assert assemblies[0].molecular_layout() is layout
    assert result.root.mol.GetNumConformers() == conformers
    assert all(a.GetIsotope() == 0 for a in layout.molecule.GetAtoms())
    assert json.loads(json.dumps(result.to_dict())) == before


def test_filtered_graphs_keep_layout_context():
    result, *_ = layered_result()
    assembly = result.assembly_graph
    expected = assembly.monomer_positions()
    filtered = assembly.drop_unassigned().filtered_by_root_bond_elements()
    assert filtered.monomer_positions() == expected
    for component in filtered.connected_components():
        assert all(position == expected[node] for node, position in component.monomer_positions().items())


def test_detached_assembly_accepts_explicit_layout_and_spring_remains_available():
    result, *_ = layered_result()
    assembly = AssemblyGraph.from_dict(result.assembly_graph.to_dict())
    with pytest.raises(ValueError, match='original root molecule'):
        assembly.monomer_positions()
    layout = MoleculeLayout.from_mol(result.root.mol)
    assert assembly.monomer_positions(layout) == result.assembly_graph.monomer_positions(layout)
    fig, ax = plt.subplots()
    try:
        assembly.draw(layout='spring', ax=ax)
    finally:
        plt.close(fig)


def test_zero_heavy_atom_depiction_does_not_invent_a_node_position():
    result = make_result('[2H][2H]', [])
    fig, ax = plt.subplots()
    try:
        result.assembly_graph.draw(ax=ax, show_unassigned=True)
        assert any(text.get_text() == 'No original heavy atoms to display' for text in ax.texts)
        assert result.assembly_graph.monomer_positions() == {}
    finally:
        plt.close(fig)


def test_daptomycin_layout_survives_loading_and_atom_renumbering(ruleset, tmp_path):
    from .data.integration_retromol import CASES

    smiles = next(smiles for name, smiles, *_ in CASES if name == 'daptomycin')
    result = run_retromol(Submission(smiles), ruleset)
    loaded = Result.from_dict(json.loads(json.dumps(result.to_dict())))
    reordered = Chem.RenumberAtoms(result.root.mol, list(reversed(range(result.root.mol.GetNumAtoms()))))
    layout = MoleculeLayout.from_mol(result.root.mol)
    for molecule in [loaded.root.mol, reordered]:
        other = MoleculeLayout.from_mol(molecule)
        assert other.bounds == pytest.approx(layout.bounds)
        assert other.atom_positions.keys() == layout.atom_positions.keys()
        for tag, position in layout.atom_positions.items():
            assert other.atom_positions[tag] == pytest.approx(position)
    assemblies = loaded.sample_assemblies(10, seed=42)
    assert len(assemblies) == 5
    positions = [assembly.monomer_positions(layout) for assembly in assemblies]
    common = set.intersection(*(set(position) for position in positions))
    assert len(common) == 13  # all thirteen amino acids stay fixed
    assert all(position[node] == positions[0][node] for position in positions for node in common)
    image = tmp_path / 'daptomycin.png'
    figures = plt.get_fignums()
    assemblies[0].draw(layout=layout, savepath=str(image), show_unassigned=True)
    assert image.read_bytes().startswith(b'\x89PNG\r\n\x1a\n')
    assert plt.get_fignums() == figures
