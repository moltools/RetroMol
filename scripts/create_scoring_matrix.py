import argparse 
import json

import numpy as np
from rdkit import Chem
from tqdm import tqdm

from retromol.retrosynthesis.chem import mol_to_fingerprint, tanimoto_similarity

parser = argparse.ArgumentParser(description='Create scoring matrix for MSA')
parser.add_argument('--monomers', type=str, help='Path to monomers file')
parser.add_argument('--output', type=str, help='Path to output csv file')
args = parser.parse_args()

monomers = json.load(open(args.monomers))
monomers = [x for x in monomers if x["name"].startswith("peptide")]
matrix = np.zeros((len(monomers), len(monomers)))

radius = 2
num_bits = 2048

for i, monomer1 in tqdm(enumerate(monomers)):
    for j, monomer2 in enumerate(monomers):
        if i == j:
            matrix[i, j] = 3 # match
        elif j < i:
            continue
        else:
            # Both are peptides
            if monomer1["name"].startswith("peptide") and monomer2["name"].startswith("peptide"):
                mol1 = Chem.MolFromSmiles(monomer1["pattern"])
                mol2 = Chem.MolFromSmiles(monomer2["pattern"])
                fp1 = mol_to_fingerprint(mol1, radius=radius, num_bits=num_bits)
                fp2 = mol_to_fingerprint(mol2, radius=radius, num_bits=num_bits)
                similarity = tanimoto_similarity(fp1, fp2)
                if similarity < 0.3:
                    score = 0
                elif similarity < 0.6:
                    score = 1
                elif similarity < 0.8:
                    score = 2
                else:
                    score = 3
                matrix[i, j] = score
                matrix[j, i] = score

labels = [x["name"].split("|")[2] for x in monomers]

# first row headers then values, save as integer values
np.savetxt(args.output, matrix, delimiter=",", header=",".join(labels), comments="", fmt='%i')
