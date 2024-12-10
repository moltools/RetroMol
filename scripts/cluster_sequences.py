#!/usr/bin/env python3
import argparse
import json
import os
import pprint
from copy import deepcopy
from multiprocessing import Pool, cpu_count, set_start_method
from typing import Callable, Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from versalign.motif import Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence

from retromol.data.motifs import MOTIFS


pretty_print = pprint.PrettyPrinter(indent=4).pprint


def smiles_to_fingerprint(smiles: str) -> np.ndarray:
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = fp_generator.GetFingerprint(mol)
    fp_arr = np.zeros((2048,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, fp_arr)

    # check if fingerprint only contains 0s and 1s
    assert np.all(np.logical_or(fp_arr == 0, fp_arr == 1))

    return fp_arr


def calc_tanimoto_similarity(arr1: np.ndarray, arr2: np.ndarray) -> float:
    # check if arrays are of same size
    assert arr1.shape == arr2.shape

    # check if both arrays only contain 0s and 1s
    assert np.all(np.logical_or(arr1 == 0, arr1 == 1))

    intersection = np.dot(arr1, arr2)
    sum_arr1 = np.sum(arr1)
    sum_arr2 = np.sum(arr2)
    score = intersection / (sum_arr1 + sum_arr2 - intersection)

    # check if score is within [0, 1]
    assert 0 <= score <= 1

    return score


parasect_substrates_as_smiles = {
    "tryptophan": r"C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N",
    "leucine": r"CC(C)C[C@@H](C(=O)O)N",
    "tyrosine": r"C1=CC(=CC=C1C[C@@H](C(=O)O)N)O",
    "serine": r"C([C@@H](C(=O)O)N)O",
    "alanine": r"C[C@@H](C(=O)O)N",
    "valine": r"CC(C)[C@@H](C(=O)O)N",
    "arginine": r"C(C[C@@H](C(=O)O)N)CN=C(N)N",
    "glutamine": r"C(CC(=O)N)[C@@H](C(=O)O)N",
    "isoleucine": r"CC[C@H](C)[C@@H](C(=O)O)N",
    "aspartic acid": r"C([C@@H](C(=O)O)N)C(=O)O",
    "cysteine": r"C([C@@H](C(=O)O)N)S",
    "proline": r"C1C[C@H](NC1)C(=O)O",
    "phenylalanine": r"C1=CC=C(C=C1)C[C@@H](C(=O)O)N",
    "asparagine": r"C([C@@H](C(=O)O)N)C(=O)N",
    "4-hydroxyphenylglycine": r"C1=CC(=CC=C1[C@@H](C(=O)O)N)O",
    "2,4-diaminobutyric acid": r"C(CN)[C@@H](C(=O)O)N",
    "beta-hydroxytyrosine": r"C1=CC(=CC=C1[C@H]([C@@H](C(=O)O)N)O)O",
    "glutamic acid": r"C(CC(=O)O)[C@@H](C(=O)O)N",
    "glycine": r"NCC(=O)O",
    "3,5-dihydroxyphenylglycine": r"N[C@H](C(=O)O)c1cc(O)cc(O)c1",
    "pipecolic acid": r"C1CCN[C@@H](C1)C(=O)O",
    "threonine": r"C[C@H]([C@@H](C(=O)O)N)O",
    "N5-formyl-N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CN(C=O)O",
    "N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CNO",
    "ornithine": r"C(C[C@@H](C(=O)O)N)CN",
    "beta-alanine": r"NCCC(=O)O",
    "histidine": r"C1=C(NC=N1)C[C@@H](C(=O)O)N",
    "D-alanine": r"C[C@H](C(=O)O)N",
    "2-aminoadipic acid": r"C(C[C@@H](C(=O)O)N)CC(=O)O",
    "lysine": r"C(CCN)C[C@@H](C(=O)O)N",
    "2-aminoisobutyric acid": r"O=C(O)C(N)(C)C",
    "2,3-dihydroxybenzoic acid": r"C1=CC(=C(C(=C1)O)O)C(=O)O",
    "salicylic acid": r"C1=CC=C(C(=C1)C(=O)O)O",
    "anthranillic acid": r"C1=CC=C(C(=C1)C(=O)O)N",
}


polyketide_motif_identities = [
    "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11",
    "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11",
    "C1", "C2", "C4",
    "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
]


def format_polyketide_tautomer_identity(name: str) -> str:
    if "-tautomer" in name:
        return name.replace("-tautomer", "")
    return name


class SequenceMotif(Motif):
    def __init__(self, name: str) -> None:
        self._name = name

    def __eq__(self, other):
        # check if other is of same type
        if not isinstance(other, SequenceMotif):
            return False

        # check if names are equal
        return self._name == other._name
    
    def __str__(self):
        return f"<{self._name}>"
    
    @property
    def name(self) -> str:
        return self._name
    


polyketide_smarts = Chem.MolFromSmarts("[OH]SC~CC(=O)[OH]")
alpha_amino_acid_smarts = Chem.MolFromSmarts("[NH2]CC(=O)[OH]")
alpha_amino_acid_smart_proline_like = Chem.MolFromSmarts("[NH1]CC(=O)[OH]")
beta_amino_acid_smarts = Chem.MolFromSmarts("[NH2]CCC(=O)[OH]")

def match_other_motif(smiles: str) -> Optional[str]:
    mol = Chem.MolFromSmiles(smiles)

    if mol.HasSubstructMatch(polyketide_smarts):
        return "polyketide"
    if mol.HasSubstructMatch(alpha_amino_acid_smarts):
        return "amino acid"
    if mol.HasSubstructMatch(alpha_amino_acid_smart_proline_like):
        return "amino acid"
    if mol.HasSubstructMatch(beta_amino_acid_smarts):
        return "amino acid"
    
    return None


def parse_motif(name: Optional[str], smiles: str, parasect_substrates_as_fingerprint: Dict[str, np.ndarray], retromol_motif_as_fingerprint: Dict[str, np.ndarray]) -> Motif:
    # unidentified motifs are named 'other'
    if name is None:
        return SequenceMotif("other")

    # check for '-tautomer' suffix in polyketide motif identities
    name = format_polyketide_tautomer_identity(name)

    # parse name
    if name in polyketide_motif_identities:
        return SequenceMotif(name)

    # check if it is substrate known by parasect, check by comparing fingerprints
    motif_fp = retromol_motif_as_fingerprint[name]
    for substrate_name, substrate_fp in parasect_substrates_as_fingerprint.items():
        score = calc_tanimoto_similarity(motif_fp, substrate_fp)
        if score == 1.0:
            return SequenceMotif(substrate_name)
        
    # check if we can identify a backbone in the motif
    if name := match_other_motif(smiles):
        return SequenceMotif(name)
    
    # unidentified motifs are named 'other'
    return SequenceMotif("other")


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="retromol results folder")
    parser.add_argument("-s", type=str, required=True, help="scoring matrix file")
    parser.add_argument("-o", type=str, required=True, help="output folder")
    parser.add_argument("-c", type=int, default=cpu_count(), help="number of cpus to use")
    return parser.parse_args()


def parse_retromol_results(results_folder: str, parasect_substrates_as_fingerprint: Dict[str, np.ndarray], retromol_motif_as_fingerprint: Dict[str, np.ndarray]) -> Tuple[List[str], List[str], List[Sequence]]:
    # instantiate data structures
    record_names = []
    record_smiles_strings = []
    record_sequences = []

    # parse over subfolders, select for put.json file in each subfolder
    for subdir in tqdm(os.listdir(results_folder), desc="parsing retromol results...", leave=True):
        data_path = os.path.join(results_folder, subdir, "out.json")

        # check if file exists, if not the item failed to parse
        if not os.path.exists(data_path):
            continue

        with open(data_path, "r") as f:
            data = json.load(f)

        # parse data
        record_name = data["name"]
        record_smiles_string = data["smiles"]
        record_coverage = data["coverage_score"]
        encoding_to_smiles = data["encoding_to_smiles"]
        encoding_to_identity = data["encoding_to_identity"]
        encoding_to_motif_codes = data["encoding_to_motif_codes"]

        # parse out sequences
        for node, motif_codes in encoding_to_motif_codes.items():
            for motif_code in motif_codes.values():
                motifs = []
                for motif in motif_code[::-1]:
                    motif_name = motif["identity"]
                    motif_smiles = motif["smiles"]
                    parsed_motif = parse_motif(motif_name, motif_smiles, parasect_substrates_as_fingerprint, retromol_motif_as_fingerprint)
                    motifs.append(parsed_motif)

                # instantiate sequence
                sequence = Sequence(record_name, motifs)

                # append to data structures
                record_names.append(record_name)
                record_smiles_strings.append(record_smiles_string)
                record_sequences.append(sequence)

    # return data structures
    return record_names, record_smiles_strings, record_sequences


def parse_scoring_matrix(path_to_file: str) -> Dict[str, Dict[str, int]]:
    scoring_matrix = {}

    # top row and first column are both headers, names might contain commas as well
    # open with pandas and parse from there to avoid issues with commas in names
    df = pd.read_csv(path_to_file, index_col=0)
    for name in df.columns:
        scoring_matrix[name] = df[name].to_dict()
    
    return scoring_matrix


def score_func(scoring_matrix: Dict[str, Dict[str, int]], a: Motif, b: Motif) -> int:
    return scoring_matrix[a.name][b.name]


def compute_similarity(args: Tuple[int, int, Sequence, Sequence, Dict[str, Dict[str, int]], int]) -> Tuple[int, int, int]:
    """Compute similarity for a given pair of indices."""
    # unpack arguments
    i, j, seq_a, seq_b, scoring_matrix, match_score = args

    seq_a = deepcopy(seq_a)
    seq_b = deepcopy(seq_b)

    def score_motif_similarity(a: Motif, b: Motif) -> float:
        return score_func(scoring_matrix, a, b)

    _, _, score = align_pairwise(
        seq_a=seq_a, seq_b=seq_b,
        score_func=score_motif_similarity,
        algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH,
        options={"gap_penalty": 3, "end_gap_penalty": 2},  # conserved sequences: (5, 3); variable sequences: (3, 1); balanced: (4, 2)
    )

    max_score = match_score * len(seq_a)
    score /= max_score

    return (i, j, score)


def pairwise_similarity(num_cpus: int, record_sequences: List[Sequence], scoring_matrix: Dict[str, Dict[str, int]], match_score: int) -> np.ndarray:
    # make sure to keep one cpu free
    num_cpus = min(num_cpus, cpu_count() - 1)
    print(f"using {num_cpus} cpus for parallel computation")

    # create an empty similarity matrix
    similarity_matrix = np.zeros((len(record_sequences), len(record_sequences)))

    # prepare tasks for multiprocessing
    tasks = []
    for i in tqdm(range(len(record_sequences)), desc="preparing tasks...", total=len(record_sequences), leave=True):
        for j in range(i, len(record_sequences)):
            seq_a = record_sequences[i]
            seq_b = record_sequences[j]

            # only add as tasks if the sequences are same length or at max 3 units different
            if abs(len(seq_a) - len(seq_b)) <= 3:
                tasks.append((i, j, seq_a, seq_b, scoring_matrix, match_score))

    print(f"number of tasks: {len(tasks)}")

    # use multiprocessing to parallelize the task
    pool = Pool(processes=num_cpus)
    try:
        results = list(tqdm(pool.imap(compute_similarity, tasks), total=len(tasks), desc="calculating pairwise similarities..."))
    finally:
        pool.close()  # prevent new tasks from being submitted
        pool.join()   # wait for all worker processes to finish

    # populate the similarity matrix with the results
    for i, j, score in results:
        max_score_i = match_score * len(record_sequences[i])
        max_score_j = match_score * len(record_sequences[j])
        similarity_matrix[i, j] = score / max_score_i
        similarity_matrix[j, i] = score / max_score_j

    return similarity_matrix


def main() -> None:
    do_test = False

    # disable rdkit logging
    RDLogger.DisableLog("rdApp.*")

    # convert parasect substrates to fingerprints
    parasect_substrates_as_fingerprint = {
        name: smiles_to_fingerprint(smiles) 
        for name, smiles 
        in tqdm(parasect_substrates_as_smiles.items(), desc="converting parasect substrates to fingerprints...", total=len(parasect_substrates_as_smiles), leave=True)
    }

    # convert retromol motifs to fingerprints
    retromol_motif_as_fingerprint = {
        motif["name"]: smiles_to_fingerprint(motif["smiles"])
        for motif in tqdm(MOTIFS, desc="converting retromol motifs to fingerprints...", total=len(MOTIFS), leave=True)
    }

    # parse command line arguments
    args = cli()

    # parse scoring matrix and use in scoring function
    scoring_matrix = parse_scoring_matrix(args.s)

    # get highest score on diagonal
    match_score = max([scoring_matrix[name][name] for name in scoring_matrix.keys()])
    print(match_score)

    # test pairwise alignment
    if do_test:
        def motif_code_to_seq(motif_code: List[str]) -> List[Motif]:
            motifs = []
            for motif_name in motif_code:
                parsed_motif = parse_motif(motif_name, "", parasect_substrates_as_fingerprint, retromol_motif_as_fingerprint)
                motifs.append(parsed_motif)
            return motifs

        seq1 = ["ethanoic acid", "C1", "C2", "B2", "B2", "D2", "C2", "B2", "C1", "B1", "B2", "B2"]
        seq2 = ["ethanoic acid", "C1", "C2", "B2", "B1", "D2", "D2", "B2", "C1", "B1", "B2", "C1", "C1"]
        # seq1 = ["B1", "A3", "B5"]
        # seq2 = ["B1", "A3", "B5"]
        seq1 = Sequence("test1", motif_code_to_seq(seq1))
        seq2 = Sequence("test2", motif_code_to_seq(seq2))

        def score_func(a: Motif, b: Motif) -> int:
            return scoring_matrix[a.name][b.name]

        aligned1, aligned2, score1 = align_pairwise(
            seq_a=seq1, seq_b=seq2,
            score_func=score_func,
            algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH,
            options={"gap_penalty": 3, "end_gap_penalty": 2},
        )
        _, _, score2 = compute_similarity((0, 1, seq1, seq2, scoring_matrix, match_score))

        print("aligned1:", aligned1)
        print("aligned2:", aligned2)
        for a, b in zip(aligned1, aligned2):
            print(a, b)
        print("scores:", score1, score2)

        exit(1)

    # parse sequences from retromol results folder
    record_names, record_smiles_strings, record_sequences = parse_retromol_results(args.i, parasect_substrates_as_fingerprint, retromol_motif_as_fingerprint)

    # filter out sequences that are too short (<6 units)
    ids_to_keep = []
    for i, sequence in enumerate(record_sequences):
        if len(sequence) >= 6:
            ids_to_keep.append(i)

    record_names = [record_names[i] for i in ids_to_keep]
    record_smiles_strings = [record_smiles_strings[i] for i in ids_to_keep]
    record_sequences = [record_sequences[i] for i in ids_to_keep]

    # randomly pick 100 from the dataset
    inds = np.random.choice(len(record_names), 1000, replace=False)
    record_names = [record_names[i] for i in inds]
    record_smiles_strings = [record_smiles_strings[i] for i in inds]
    record_sequences = [record_sequences[i] for i in inds]
    
    print("number of records:", len(record_names))
    print("number of smiles strings:", len(record_smiles_strings))
    print("number of sequences:", len(record_sequences))

    # calculate pairwise similarity
    similarity_matrix = pairwise_similarity(args.c, record_sequences, scoring_matrix, match_score)

    # save similarity matrix to file, as well as the record names and smiles strings
    out_tuple = (record_names, record_smiles_strings, similarity_matrix)
    out_path = os.path.join(args.o, "similarity_matrix.pkl")
    with open(out_path, "wb") as f:
        pickle.dump(out_tuple, f)


if __name__ == "__main__":
    main()
