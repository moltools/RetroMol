#!/usr/bin/env python3

import argparse

import numpy as np
import yaml

from retromol.chem.mol import smiles_to_mol
from retromol.chem.fingerprint import mol_to_morgan_fingerprint, calculate_tanimoto_similarity


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mxn", type=str, required=True, help="Path to matching rules file.")
    parser.add_argument("--out", type=str, required=True, help="Path to to output csv file.")
    return parser.parse_args()


def main() -> None:
    args = cli()

    records: list[tuple[str, str]] = []
    with open(args.mxn, "r") as fo:
        mxn = yaml.safe_load(fo)
    for x in mxn:
        name = x["name"]
        smiles = x["smiles"]
        records.append((name, smiles))
    print(len(records))

    labels = [x[0] for x in records]
    mols = [smiles_to_mol(x[1]) for x in records]
    print(len(labels), len(mols))

    fingerprints = [mol_to_morgan_fingerprint(mol, radius=2, num_bits=2048, use_chirality=False) for mol in mols]

    similarity_matrix = np.zeros((len(fingerprints), len(fingerprints)))
    for i in range(len(fingerprints)):
        for j in range(i, len(fingerprints)):
            if i == j:
                similarity_matrix[i, j] = 1.0
                continue
            elif i > j:
                continue
            similarity = calculate_tanimoto_similarity(fingerprints[i], fingerprints[j])
            similarity_matrix[i, j] = similarity
            similarity_matrix[j, i] = similarity
    print(similarity_matrix.shape)

    # round all floats in similarity matrix to 2 after decimal
    similarity_matrix = np.round(similarity_matrix, 2)

    # add "..." to labels
    labels = [f'\"{label}\"' for label in labels]

    # Write similaritym matrix to csv file with labels as column and row headers; intersection column/row top left should be 'identity'
    with open(args.out, "w") as fo:
        fo.write("identity," + ",".join(labels) + "\n")
        for i in range(len(labels)):
            fo.write(labels[i] + "," + ",".join([str(similarity_matrix[i, j]) for j in range(len(labels))]) + "\n")


if __name__ == "__main__":
    main()
