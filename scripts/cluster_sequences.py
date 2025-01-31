#!/usr/bin/env python3
import argparse
import json
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import cpu_count
from typing import Dict, List, Optional, Tuple
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
    "A", "B", "C", "D",
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
        return f"{self._name}"
    
    def __hash__(self):
        return hash(self._name)
    
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

def parse_motif(name: Optional[str], smiles: str, motif_to_fingerprint, substrate_to_fingerprint) -> Motif:
    # unidentified motifs are named 'other'
    if name is None or name == "other":
        return SequenceMotif("other")

    # check if start or end
    if name == "start" or name == "end":
        return SequenceMotif(name)

    # check for '-tautomer' suffix in polyketide motif identities
    name = format_polyketide_tautomer_identity(name)

    # parse name
    if name in polyketide_motif_identities:
        return SequenceMotif(name)

    # check if it is substrate known by parasect, check by comparing fingerprints
    motif_fp = motif_to_fingerprint[name]
    for substrate_name, substrate_fp in substrate_to_fingerprint.items():
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
    parser.add_argument("--parsed_data", type=str, default=None, help="path to parsed data file")

    # flag to say deorphanize
    parser.add_argument("--deorphanize", action="store_true", help="deorphanize motifs")
    parser.add_argument("--mibig-mapping", type=str, default=None, required=False, help="path to mibig mapping file")

    return parser.parse_args()
    
def extract_kmers(word_list, k):
    """Extract k-mers from a list of words."""
    kmers = []
    for i in range(len(word_list) - k + 1):
        kmers.append(tuple(word_list[i:i + k]))
    return kmers

def parse_retromol_results(results_folder: str, motif_to_fingerprint, substrate_to_fingerprint) -> Tuple[List[str], List[str], List[Sequence]]:
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
                motifs = ["start"]
                for motif in motif_code[::-1]:
                    motif_name = motif["identity"]
                    motif_smiles = motif["smiles"]
                    parsed_motif_string = parse_motif(motif_name, motif_smiles, motif_to_fingerprint, substrate_to_fingerprint)
                    motifs.append(parsed_motif_string)
                motifs.append("end")

                # instantiate sequence
                motifs = [SequenceMotif(motif) for motif in motifs]
                sequence = Sequence(record_name, motifs)

                # append to data structures
                record_names.append(record_name)
                record_smiles_strings.append(record_smiles_string)
                record_sequences.append(sequence)

    # return data structures
    return record_names, record_smiles_strings, record_sequences

def parse_scoring_matrix(path_to_file: str) -> Dict[Tuple[str, str], int]:
    scoring_matrix = {}

    # top row and first column are both headers, names might contain commas as well
    # open with pandas and parse from there to avoid issues with commas in names
    df = pd.read_csv(path_to_file, index_col=0)
    for name in df.columns:
        for other_name, score in df[name].items():
            scoring_matrix[(name, other_name)] = score
    
    return scoring_matrix

def score_func(scoring_matrix: Dict[str, Dict[str, int]], a: Motif, b: Motif) -> int:
    return scoring_matrix[(a.name, b.name)]

# Function to compute similarity and alignment length for a pair of indices
def compute_similarity(args):
    i, j, record_sequences, scoring_matrix, match_score = args

    seq_a = record_sequences[i]
    seq_b = record_sequences[j]

    def score_func(a, b):
        return scoring_matrix[(str(a.name), str(b.name))]

    aligned_a, aligned_b, score = align_pairwise(
        seq_a=seq_a, seq_b=seq_b,
        score_func=score_func,
        algorithm=PairwiseAlignment.SMITH_WATERMAN,
        options={"gap_penalty": 3}
    )

    alignment_length = len(aligned_a)
    normalized_score = score / (alignment_length * match_score)
    scaled_normalized_score = int(normalized_score * 100)
    assert 0 <= scaled_normalized_score <= 100, f"score: {normalized_score}, alignment_length: {alignment_length}, score: {score}, match_score: {match_score}"

    return i, j, scaled_normalized_score, alignment_length

# Parallelized similarity matrix computation with chunking
def parallel_similarity_matrix(record_sequences, scoring_matrix, match_score, num_workers=1, chunk_size=10000):
    n = len(record_sequences)
    combined_matrix = np.zeros((n, n), dtype=np.int8)

    # Determine the number of workers
    all_cpus = cpu_count()
    num_workers = min(num_workers, all_cpus - 1)
    print(f"Assigned {num_workers} workers out of {all_cpus} available.")

    # Function to generate index pairs for chunks
    def generate_index_pairs():
        for i in range(n):
            for j in range(i, n):
                yield i, j, record_sequences, scoring_matrix, match_score

    # Create a generator for index pairs
    index_pair_generator = generate_index_pairs()

    # Process in chunks
    total_pairs = n * (n + 1) // 2
    print(f"Total pairs: {total_pairs}; processing in chunks of {chunk_size} pairs.")
    completed_pairs = 0
    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        with tqdm(total=total_pairs, desc="Calculating pairwise similarity...") as pbar:
            while completed_pairs < total_pairs:
                # Create a batch of tasks
                batch_pairs = [next(index_pair_generator) for _ in range(min(chunk_size, total_pairs - completed_pairs))]

                # Submit tasks to the executor
                future_to_indices = {executor.submit(compute_similarity, args): args[:2] for args in batch_pairs}

                # Process results as they are completed
                for future in as_completed(future_to_indices):
                    i, j = future_to_indices[future]
                    try:
                        i, j, scaled_score, alignment_length = future.result()
                        # Update the combined matrix
                        combined_matrix[i, j] = scaled_score  # Upper triangle
                        combined_matrix[j, i] = alignment_length  # Lower triangle
                    except Exception as e:
                        print(f"Error processing pair ({i}, {j}): {e}")
                    pbar.update(1)

                completed_pairs += len(batch_pairs)

    return combined_matrix


# Function to compute similarity for a row
def compute_similarity_row(args):
    i, record_sequences, scoring_matrix, match_score = args
    n = len(record_sequences)
    row_scores = np.zeros(n, dtype=np.int8)  # Store scaled scores
    row_lengths = np.zeros(n, dtype=np.int8)  # Store alignment lengths

    seq_a = record_sequences[i]

    def score_func(a, b):
        return scoring_matrix[(str(a.name), str(b.name))]

    for j in range(i, n):  # Compute only upper triangle
        seq_b = record_sequences[j]

        aligned_a, aligned_b, score = align_pairwise(
            seq_a=seq_a, seq_b=seq_b,
            score_func=score_func,
            algorithm=PairwiseAlignment.SMITH_WATERMAN,
            options={"gap_penalty": 3}
        )

        alignment_length = len(aligned_a)
        normalized_score = score / (alignment_length * match_score)
        scaled_normalized_score = int(normalized_score * 100)
        assert 0 <= scaled_normalized_score <= 100, (
            f"score: {normalized_score}, alignment_length: {alignment_length}, score: {score}, match_score: {match_score}"
        )

        # Store results for the current row
        row_scores[j] = scaled_normalized_score
        row_lengths[j] = alignment_length

    return i, row_scores, row_lengths


# Parallelized similarity matrix computation
def parallel_similarity_matrix_by_row(record_sequences, scoring_matrix, match_score, num_workers=1):
    n = len(record_sequences)
    combined_matrix = np.zeros((n, n), dtype=np.int8)  # Matrix to store scaled scores and lengths

    # Determine the number of workers
    all_cpus = cpu_count()
    num_workers = min(num_workers, all_cpus - 1)
    print(f"Assigned {num_workers} workers out of {all_cpus} available.")

    # Generate arguments for each row
    row_args = [i for i in range(n)]

    with ProcessPoolExecutor(max_workers=num_workers) as executor:
        with tqdm(total=n, desc="Calculating pairwise similarity by row...") as pbar:
            futures = {executor.submit(compute_similarity_row, (i, record_sequences, scoring_matrix, match_score)): i for i in row_args}

            for future in as_completed(futures):
                i = futures[future]
                try:
                    i, row_scores, row_lengths = future.result()

                    # Update the combined matrix
                    combined_matrix[i, i:] = row_scores[i:]  # Upper triangle
                    combined_matrix[i:, i] = row_lengths[i:]  # Lower triangle
                except Exception as e:
                    print(f"Error processing row {i}: {e}")
                pbar.update(1)

    return combined_matrix

def main() -> None:
    do_test = False
    retrieval = False
    similarity = True

    # disable rdkit logging
    RDLogger.DisableLog("rdApp.*")

    # parse command line arguments
    args = cli()

    # parse scoring matrix and use in scoring function
    scoring_matrix = parse_scoring_matrix(args.s)

    # get highest score
    match_score = max(scoring_matrix.values())

    # convert parasect substrates to fingerprints
    substrate_to_fingerprint = {
        name: smiles_to_fingerprint(smiles) 
        for name, smiles 
        in tqdm(parasect_substrates_as_smiles.items(), desc="converting parasect substrates to fingerprints...", total=len(parasect_substrates_as_smiles), leave=True)
    }

    # convert retromol motifs to fingerprints
    motif_to_fingerprint = {
        motif["name"]: smiles_to_fingerprint(motif["smiles"])
        for motif in tqdm(MOTIFS, desc="converting retromol motifs to fingerprints...", total=len(MOTIFS), leave=True)
    }

    # parse sequences from retromol results folder
    if not args.parsed_data:
        record_names, record_smiles_strings, record_sequences = parse_retromol_results(args.i, motif_to_fingerprint, substrate_to_fingerprint)
        # pickle data
        out_tuple = (record_names, record_smiles_strings, record_sequences)
        out_path = os.path.join(args.o, "parsed_data.pkl")
        with open(out_path, "wb") as f:
            pickle.dump(out_tuple, f)
    else:
        # load from pickle file
        with open(args.parsed_data, "rb") as f:
            record_names, record_smiles_strings, record_sequences = pickle.load(f)

    # filter out sequences that are too short (<6 units)
    ids_to_keep = []
    for i, sequence in enumerate(record_sequences):
        if len(sequence) >= 6:
            ids_to_keep.append(i)

    record_names = [record_names[i] for i in ids_to_keep]
    record_smiles_strings = [record_smiles_strings[i] for i in ids_to_keep]
    record_sequences = [record_sequences[i] for i in ids_to_keep]
    
    print("number of records:", len(record_names))
    print("number of smiles strings:", len(record_smiles_strings))
    print("number of sequences:", len(record_sequences))

    if args.deorphanize:
        # parse mibig mapping
        mibig_mapping = {}
        with open(args.mibig_mapping, "r") as f:
            f.readline()
            for line in f:
                parts = line.strip().split(",")
                mibig_accession = parts[1]
                npaids = parts[2].split("|")
                mibig_mapping[mibig_accession] = npaids

        mibig_accession = "BGC0000055"
        target_npaids = mibig_mapping[mibig_accession]
        print("target npaids:", target_npaids)
        query_string = ["start", "other", "B", "B", "B", "D", "B", "B", "end"]
        query_motifs = [parse_motif(name, "", motif_to_fingerprint, substrate_to_fingerprint) for name in query_string]

        def score_func(a: SequenceMotif, b: SequenceMotif) -> int:
            return scoring_matrix[(str(a.name), str(b.name))]

        # find best 10 best matches with pairwise similarity
        query_seq = Sequence("query", query_motifs)
        query_scores = []
        for i, seq in tqdm(enumerate(record_sequences)):
            aligned_a, aligned_b, score_a = align_pairwise(
                seq_a=query_seq, seq_b=seq,
                score_func=score_func,
                algorithm=PairwiseAlignment.SMITH_WATERMAN,
                options={"gap_penalty": 2}
            )
            query_scores.append((i, score_a, seq, aligned_a, aligned_b))
        
        # sort on the two scores, makes sure both scores are as high as possible
        query_scores.sort(key=lambda x: x[1], reverse=True)

        # print position of highest ranking target npaid
        for j, (i, score_a, seq, aligned_a, aligned_b) in enumerate(query_scores):
            if seq._identifier in target_npaids:
                print(seq._identifier, target_npaids)
                print("target npaid found at position:", j, "out of", len(record_sequences))
                # break
        
        print("best hits:")
        for j, (i, score_a, seq, aligned_a, aligned_b) in enumerate(query_scores[:10]):
            print(i, record_names[i], score_a, "|".join([str(m) for m in aligned_a._motifs]), "|".join([str(m) for m in aligned_b._motifs]))

        exit(1)

    # test pairwise alignment
    if do_test:
        def motif_code_to_seq(motif_code: List[str]) -> List[Motif]:
            motifs = []
            for motif_name in motif_code:
                parsed_motif = parse_motif(motif_name, "", motif_to_fingerprint, substrate_to_fingerprint)
                motifs.append(parsed_motif)
            return motifs

        seq1 = ["ethanoic acid", "C1", "C2", "B2", "B1", "D2", "D2", "B2", "C1", "B1", "B2", "C1", "C1"]
        # seq2 = ["ethanoic acid", "C1", "C2", "B2", "B2", "D2", "C2", "B2", "C1", "B1", "B2", "B2"]
        # seq1 = ["B6", "B2", "A2", "D6", "B2", "B2"]
        seq2 = ["B2", "B2", "A2", "D2", "B2", "B2"]
        seq1 = Sequence("test1", motif_code_to_seq(seq1))
        seq2 = Sequence("test2", motif_code_to_seq(seq2))

        def score_func(a: Motif, b: Motif) -> int:
            return scoring_matrix[(a.name, b.name)]

        # aligned1, aligned2, score = align_pairwise(
        #     seq_a=seq1, seq_b=seq2,
        #     score_func=score_func,
        #     algorithm=PairwiseAlignment.NEEDLEMAN_WUNSCH,
        #     options={"gap_penalty": 3, "end_gap_penalty": 2},
        # )

        aligned1, aligned2, score = align_pairwise(
            seq_a=seq1, seq_b=seq2,
            score_func=score_func,
            algorithm=PairwiseAlignment.SMITH_WATERMAN,
            options={"gap_penalty": 3},
        )

        print("aligned1:", aligned1)
        print("aligned2:", aligned2)
        for a, b in zip(aligned1, aligned2):
            print(a, b)
        # max_score = match_score * len(seq1)
        max_score = match_score * len(aligned1)
        print("scores:", score, score / max_score)

        exit(1)

    if retrieval:
        # query_string = ["start", "B", "B", "B", "D", "B", "B", "end"]
        # query_string = ["start", "tryptophan", "asparagine", "aspartic acid", "threonine", "glycine", "ornithine", "aspartic acid", "alanine", "aspartic acid", "glycine", "serine", "3-methylglutamic acid", "kynurenine", "end"]
        query_string = ["C4", "C1", "B1", "C1", "end"] # anguinomycin/ratjadon(e)/leptofuranin/leptomycin/leptolstatin/roimatacene(NPA014832)/Nafuredin B(NPA022519) nuclear charge
        query_motifs = [parse_motif(name, "", motif_to_fingerprint, substrate_to_fingerprint) for name in query_string]

        def score_func(a: SequenceMotif, b: SequenceMotif) -> int:
            return scoring_matrix[(str(a.name), str(b.name))]

        # find best 10 best matches with pairwise similarity
        query_seq = Sequence("query", query_motifs)
        query_scores = []
        for i, seq in tqdm(enumerate(record_sequences)):
            aligned_a, aligned_b, score_a = align_pairwise(
                seq_a=query_seq, seq_b=seq,
                score_func=score_func,
                algorithm=PairwiseAlignment.SMITH_WATERMAN,
                options={"gap_penalty": 3}
            )
            query_scores.append((i, score_a, seq, aligned_a, aligned_b))
        
        # sort on the two scores, makes sure both scores are as high as possible
        query_scores.sort(key=lambda x: x[1], reverse=True)
        
        print("best hits:")
        for i, score_a, seq, aligned_a, aligned_b in query_scores[:20]:
            print(record_names[i], score_a, "|".join([str(m) for m in aligned_a._motifs]), "|".join([str(m) for m in aligned_b._motifs]))

        exit(1)

    if similarity:
        
        print(f"Calculating similarity matrix with {args.c} workers...")
        # combined_matrix = parallel_similarity_matrix(record_sequences, scoring_matrix, match_score, num_workers=args.c)
        combined_matrix = parallel_similarity_matrix_by_row(record_sequences, scoring_matrix, match_score, num_workers=args.c)

        # # calculate pairwise similarity
        # combined_matrix = np.zeros((len(record_sequences), len(record_sequences)), dtype=np.int8)
        # for i in tqdm(range(len(record_sequences)), desc="calculating pairwise similarity...", leave=True):
        #     for j in tqdm(range(i, len(record_sequences)), desc="calculating pairwise similarity...", leave=False):
        #         seq_a = record_sequences[i]
        #         seq_b = record_sequences[j]

        #         def score_func(a: SequenceMotif, b: SequenceMotif) -> int:
        #             return scoring_matrix[(str(a.name), str(b.name))]

        #         aligned_a, aligned_b, score = align_pairwise(
        #             seq_a=seq_a, seq_b=seq_b,
        #             score_func=score_func,
        #             algorithm=PairwiseAlignment.SMITH_WATERMAN,
        #             options={"gap_penalty": 3}
        #         )

        #         # calculate similarity score
        #         alignment_length = len(aligned_a)
        #         normalized_score = score / (alignment_length * match_score)
        #         scaled_normalized_score = int(normalized_score * 100)
        #         assert 0 <= scaled_normalized_score <= 100, f"score: {normalized_score}, alignment_length: {alignment_length}, score: {score}, match_score: {match_score}"

        #         # score values score and alignment length in same combined matrix, upper and lower triangles
        #         combined_matrix[i, j] = scaled_normalized_score # triangle: upper
        #         combined_matrix[j, i] = alignment_length # triangle: lower

        out_tuple = (record_names, record_smiles_strings, combined_matrix)
        out_path = os.path.join(args.o, "similarity_matrix.pkl")
        with open(out_path, "wb") as f:
            pickle.dump(out_tuple, f)

if __name__ == "__main__":
    main()
