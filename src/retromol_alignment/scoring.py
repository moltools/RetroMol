"""Scoring functionalities for sequence alignment."""

from typing import Sequence

import numpy as np

from retromol.chem.mol import smiles_to_mol
from retromol.chem.fingerprint import mol_to_morgan_fingerprint, calculate_tanimoto_similarity
from retromol_alignment.aligner import substitution_matrices


def create_tanimoto_scoring_matrix(
    records: Sequence[tuple[str, str]],
    radius: int = 2,
    num_bits: int = 2048,
    stereochemistry: bool = False,
    self_score_tokens: Sequence[str] | None = None,
    self_score: float = 1.0,
) -> substitution_matrices.Array:
    names: list[str] = []
    fps = []

    for name, smiles in records:
        mol = smiles_to_mol(smiles)
        fp = mol_to_morgan_fingerprint(mol, radius=radius, num_bits=num_bits, use_chirality=stereochemistry)
        names.append(name)
        fps.append(fp)

    if self_score_tokens is not None:
        names.extend(self_score_tokens)

    n_mols = len(fps)
    n = len(names)

    data = np.zeros((n, n), dtype=np.float64)

    for i in range(n_mols):
        for j in range(i, n_mols):
            if i == j:
                score = 1.0
            else:
                score = calculate_tanimoto_similarity(fps[i], fps[j])

            data[i, j] = score
            data[j, i] = score

    for i in range(n_mols, n):
        data[i, i] = self_score

    alphabet = tuple(names)

    return substitution_matrices.Array(alphabet, 2, data, np.float64)
