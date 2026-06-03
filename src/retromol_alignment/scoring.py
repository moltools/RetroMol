"""Scoring functionalities for sequence alignment."""

from typing import Sequence

import numpy as np

from retromol.chem.mol import smiles_to_mol
from retromol.chem.fingerprint import mol_to_morgan_fingerprint, calculate_tanimoto_similarity
from retromol_alignment.aligner import substitution_matrices


PK_A_SUBTYPES = [
    "A1",
    "A2", "A2^R", "A2^S",
    "A3",
    "A4", "A4^R", "A4^S",
    "A5", "A5^R", "A5^S",
    "A6", "A6^R", "A6^S",
    "A7", "A7^R", "A7^S",
    "A8",
    "A9", "A9^R", "A9^S",
    "A10", "A10^R", "A10^S",
    "A11",
]
PK_B_SUBTYPES = [
    "B1", "B^R1", "B^S1",
    "B2", "B^R2^R", "B^R2^S", "B^S2^R", "B^S2^S", "B^R2", "B^S2",
    "B3", "B^R3", "B^S3",
    "B4", "B^R4^R", "B^R4^S", "B^S4^R", "B^S4^S",
    "B5", "B^R5^R", "B^R5^S", "B^S5^R", "B^S5^S",
    "B6", "B^R6^R", "B^R6^S", "B^S6^R", "B^S6^S",
    "B7", "B^R7^R", "B^R7^S", "B^S7^R", "B^S7^S",
    "B8", "B^R8", "B^S8",
    "B9", "B^R9^R", "B^R9^S", "B^S9^R", "B^S9^S",
    "B10", "B^R10^R", "B^R10^S", "B^S10^R", "B^S10^S",
    "B11",
    "B12", "B^R12^R", "B^R12^S", "B^S12^R", "B^S12^S",
    "B13", "B^R13", "B^S13",
]
PK_C_SUBTYPES = [
    "C1",
    "C2",
    "C4",
    "C7",
    "C13",
]
PK_D_SUBTYPES = [
    "D1",
    "D2", "D2^R", "D2^S",
    "D3",
    "D4", "D4^R", "D4^S",
    "D5", "D5^R", "D5^S",
    "D6", "D6^R", "D6^S",
    "D7", "D7^R", "D7^S",
    "D8", "D^R8", "D^S8",
    "D10", "D^R10^R", "D^R10^S", "D^S10^R", "D^S10^S",
    "D11",
    "D14",
    "D15",
    "D16", "D16^R", "D16^S",
    "D17", "D^R17", "D^S17",
]
HARDCODED_PK_SCORING = []
HARDCODED_PK_SCORING.extend([("A", pk_subtype, 1.0) for pk_subtype in PK_A_SUBTYPES])
HARDCODED_PK_SCORING.extend([("B", pk_subtype, 1.0) for pk_subtype in PK_B_SUBTYPES])
HARDCODED_PK_SCORING.extend([("C", pk_subtype, 1.0) for pk_subtype in PK_C_SUBTYPES])
HARDCODED_PK_SCORING.extend([("D", pk_subtype, 1.0) for pk_subtype in PK_D_SUBTYPES])



def create_tanimoto_scoring_matrix(
    records: Sequence[tuple[str, str]],
    radius: int = 2,
    num_bits: int = 2048,
    stereochemistry: bool = False,
    self_score_tokens: Sequence[str] | None = None,
    self_score: float = 1.0,
    hardcoded_scores: list[tuple[str, str, float]] | None = None,
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

    # Finally, supplement scores with hardcoded scores, if any
    for a, b, score in hardcoded_scores or []:
        if a not in names or b not in names:
            raise ValueError(f"Hardcoded score tokens must be in the alphabet: {a}, {b}")
        i = names.index(a)
        j = names.index(b)
        data[i, j] = score
        data[j, i] = score

    alphabet = tuple(names)

    return substitution_matrices.Array(alphabet, 2, data, np.float64)
