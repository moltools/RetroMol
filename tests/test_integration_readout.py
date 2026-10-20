"""Exact, ordered sequence expectations for cyclic peptide/PK/hybrid inputs."""

import json

import pytest
from rdkit import Chem

from retromol.model.assembly_graph import AssemblyGraph
from retromol.model.readout import assembly_to_primary_sequence, split_named_paths
from retromol.model.result import Result

from .data.integration_retromol import CASES as ASSEMBLY_CASES
from .helpers_retromol import parse_compound

EPOTHILONE = (
    r'C[C@H]1CCC[C@@]2(C)O[C@H]2C[C@H](OC(=O)C[C@H](O)C(C)(C)'
    r'C(=O)[C@H](C)[C@H]1O)C(\C)=C\C1=CSC(C)=N1'
)
ERYTHROMYCIN_C = (
    'CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]('
    '[C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)'
    'O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O'
)
NEOPELTOLIDE = (
    r'O=C(/C=C\CCC1=COC(/C=C\CNC(OC)=O)=N1)O[C@@H]2C[C@@H]('
    r'C[C@@H](C)C[C@H](OC)C[C@H](CCC)OC(C3)=O)O[C@@H]3C2'
)
THIOCORALINE = (
    'CSCC1N(C)C(=O)C2CSSCC(N(C)C(=O)CNC(=O)C(CSC1=O)NC(=O)'
    'C1=NC3=CC=CC=C3C=C1O)C(=O)N(C)C(CSC)C(=O)SCC(NC(=O)'
    'C1=C(O)C=C3C=CC=CC3=N1)C(=O)NCC(=O)N2C'
)


READOUT_CASES = [
    pytest.param(
        dict(ASSEMBLY_CASES)['actinoquinoline B'],
        [['3-hydroxyquinaldic acid', '5-hydroxylysine'], ['acetic acid starter', 'D8']],
        [['3-hydroxyquinaldic acid', '5-hydroxylysine'], ['acetic acid starter', 'D8']],
        id='actinoquinoline-B',
    ),
    pytest.param(
        dict(ASSEMBLY_CASES)['chaetoglobosin A'],
        [['acetic acid starter', 'C2', 'C1', 'C1', 'D2', 'C2', 'B11', 'C1', 'A1', 'tryptophan'],
         ['oxidation']],
        [['acetic acid starter', 'C2', 'C1', 'C^E1', 'D2^R', 'C^E2', 'B11', 'C^E1', 'A1', 'tryptophan^L'],
         ['oxidation']],
        id='chaetoglobosin-A',
    ),
    pytest.param(
        dict(ASSEMBLY_CASES)['chlorothricin'],
        [['acetic acid starter', 'C13', 'C1', 'C1', 'D1', 'D1', 'C1', 'C1', 'B1', 'D1', 'C1'],
         ['acetic acid starter', 'C1', 'C1', 'A1'],
         ['glycosylation', 'glycosylation'], ['chlorination'], ['glyceric acid'],
         ['glycolic acid'], ['methylation']],
        [['acetic acid starter', 'C13', 'C1', 'C1', 'D1', 'D1', 'C1', 'C1', 'B^R1', 'D1', 'C1'],
         ['acetic acid starter', 'C1', 'C1', 'A1'],
         ['glycosylation', 'glycosylation'], ['chlorination'], ['glyceric acid'],
         ['glycolic acid'], ['methylation']],
        id='chlorothricin',
    ),
    pytest.param(
        dict(ASSEMBLY_CASES)['cremimycin'],
        [['propanoic acid starter', 'D1', 'D1', 'A', 'A1', 'D1', 'B1', 'A1', 'B1', 'C1', 'C2', 'C1', 'C1'],
         ['methylation', 'glycosylation', 'oxidation']],
        [['propanoic acid starter', 'D1', 'D1', 'A', 'A1', 'D1', 'B1', 'A1', 'B1', 'C1', 'C^E2', 'C^Z1', 'C^Z1'],
         ['methylation', 'glycosylation', 'oxidation']],
        id='cremimycin',
    ),
    pytest.param(
        dict(ASSEMBLY_CASES)['enterobactin'],
        [['2,3-dihydroxybenzoic acid', 'serine'],
         ['2,3-dihydroxybenzoic acid', 'serine'],
         ['2,3-dihydroxybenzoic acid', 'serine']],
        [['2,3-dihydroxybenzoic acid', 'serine'],
         ['2,3-dihydroxybenzoic acid', 'serine'],
         ['2,3-dihydroxybenzoic acid', 'serine']],
        id='enterobactin',
    ),
    pytest.param(
        THIOCORALINE,
        [['3-hydroxyquinaldic acid', 'cysteine', 'glycine', 'N-methylcysteine', 'N-methylcysteine'],
         ['3-hydroxyquinaldic acid', 'cysteine', 'glycine', 'N-methylcysteine', 'N-methylcysteine'],
         ['methylation'], ['methylation']],
        [['3-hydroxyquinaldic acid', 'cysteine', 'glycine', 'N-methylcysteine', 'N-methylcysteine'],
         ['3-hydroxyquinaldic acid', 'cysteine', 'glycine', 'N-methylcysteine', 'N-methylcysteine'],
         ['methylation'], ['methylation']],
        id='thiocoraline',
    ),
    pytest.param(
        NEOPELTOLIDE,
        [['acetic acid starter', 'D1', 'B1', 'B1', 'D8', 'B1', 'B1', 'C1'],
         ['glycine', 'C1', 'serine', 'D1', 'C1'],
         ['carbonic acid', 'methylation'], ['methylation']],
        [['acetic acid starter', 'D1', 'B^S1', 'B^R1', 'D^R8', 'B^R1', 'B^R1', 'C1'],
         ['glycine', 'C^Z1', 'serine', 'D1', 'C^Z1'],
         ['carbonic acid', 'methylation'], ['methylation']],
        id='neopeltolide',
    ),
    pytest.param(
        EPOTHILONE,
        [['acetic acid starter', 'cysteine', 'C2', 'B1', 'C2', 'D1', 'D2', 'B2', 'A3', 'B1'],
         ['oxidation']],
        [['acetic acid starter', 'cysteine', 'C^E2', 'B^R1', 'C2', 'D1', 'D2^R', 'B^S2^S', 'A3', 'B^R1'],
         ['oxidation']],
        id='epothilone',
    ),
    pytest.param(
        dict(ASSEMBLY_CASES)['daptomycin'],
        [['acetic acid starter', 'D1', 'D1', 'D1', 'D1', 'tryptophan', 'asparagine',
          'aspartic acid', 'threonine', 'glycine', 'ornithine', 'aspartic acid', 'alanine',
          'aspartic acid', 'glycine', 'serine', '3-methylglutamic acid', 'kynurenine']],
        [['acetic acid starter', 'D1', 'D1', 'D1', 'D1', 'tryptophan^L', 'asparagine^D',
          'aspartic acid^L', 'threonine', 'glycine', 'ornithine^L', 'aspartic acid^L', 'alanine^D',
          'aspartic acid^L', 'glycine', 'serine^D', '3-methylglutamic acid', 'kynurenine']],
        id='daptomycin',
    ),
    pytest.param(
        ERYTHROMYCIN_C,
        [['propanoic acid starter', 'B6', 'B2', 'A2', 'D6', 'B2', 'B2'],
         ['glycosylation'], ['glycosylation']],
        [['propanoic acid starter', 'B^R6^S', 'B^S2^S', 'A2^R', 'D6^S', 'B^S2^S', 'B^S2^S'],
         ['glycosylation'], ['glycosylation']],
        id='erythromycin-C',
    ),
]


@pytest.mark.parametrize('smiles, expected_plain, expected_stereo', READOUT_CASES)
@pytest.mark.parametrize('stereo', [False, True], ids=['plain', 'stereo'])
@pytest.mark.parametrize('reverse_atoms', [False, True], ids=['input-order', 'reversed-atoms'])
def test_integration_linear_readout(smiles, expected_plain, expected_stereo, stereo, reverse_atoms,
                                    ruleset, stereo_ruleset):
    if reverse_atoms:
        mol = Chem.MolFromSmiles(smiles)
        mol = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
        smiles = Chem.MolToSmiles(mol, canonical=False)
    result = parse_compound(smiles, stereo_ruleset if stereo else ruleset)
    expected = expected_stereo if stereo else expected_plain
    assembly_before = result.assembly_graph.to_dict()

    readout = result.linear_readout
    assert readout.primary_sequence_parts() == expected
    assert split_named_paths(readout.primary_sequence(), key=lambda name: name) == expected
    assert result.calculate_coverage() == 1.0
    assert result.assembly_graph.to_dict() == assembly_before

    # Saved results and the assembly-only API used by sampled CLI exports must
    # preserve direction and roles without rematching or access to a live root.
    restored = Result.from_dict(json.loads(json.dumps(result.to_dict())))
    assert restored.linear_readout.primary_sequence_parts() == expected
    assembly = AssemblyGraph.from_dict(json.loads(json.dumps(assembly_before)))
    assert assembly_to_primary_sequence(assembly) == expected

    path_nodes = [node.enc for path in readout.paths for node in path]
    identified = {node.enc for node in result.assembly_graph.monomer_nodes() if node.is_identified}
    assert len(path_nodes) == len(set(path_nodes))
    assert set(path_nodes) == identified
