#!/usr/bin/env python3
import argparse 
import pandas as pd
import numpy as np
from tqdm import tqdm
from cluster_sequences import SequenceMotif
from versalign.motif import Motif, Gap
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence
from versalign.msa import multiple_sequence_alignment
from matplotlib import pyplot as plt

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cluster-assignments", type=str, required=True, help="path to cluster assignments file")
    parser.add_argument("--scoring-matrix", type=str, required=True, help="path to scoring matrix")
    parser.add_argument("--output", type=str, required=True, help="path to output folder")
    return parser.parse_args()

def parse_scoring_matrix(path_to_file: str) -> dict:
    scoring_matrix = {}
    df = pd.read_csv(path_to_file, index_col=0)
    for name in df.columns:
        for other_name, score in df[name].items():
            scoring_matrix[(name, other_name)] = score
    return scoring_matrix

def visualize_alignment(sequences: list, output: str) -> None:

    # Define a mapping of residues to colors
    colors = {
        # "tryptophan": r"C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N",
        # "leucine": r"CC(C)C[C@@H](C(=O)O)N",
        # "tyrosine": r"C1=CC(=CC=C1C[C@@H](C(=O)O)N)O",
        # "serine": r"C([C@@H](C(=O)O)N)O",
        # "alanine": r"C[C@@H](C(=O)O)N",
        # "valine": r"CC(C)[C@@H](C(=O)O)N",
        # "arginine": r"C(C[C@@H](C(=O)O)N)CN=C(N)N",
        # "glutamine": r"C(CC(=O)N)[C@@H](C(=O)O)N",
        # "isoleucine": r"CC[C@H](C)[C@@H](C(=O)O)N",
        # "aspartic acid": r"C([C@@H](C(=O)O)N)C(=O)O",
        # "cysteine": r"C([C@@H](C(=O)O)N)S",
        # "proline": r"C1C[C@H](NC1)C(=O)O",
        # "phenylalanine": r"C1=CC=C(C=C1)C[C@@H](C(=O)O)N",
        # "asparagine": r"C([C@@H](C(=O)O)N)C(=O)N",
        # "4-hydroxyphenylglycine": r"C1=CC(=CC=C1[C@@H](C(=O)O)N)O",
        # "2,4-diaminobutyric acid": r"C(CN)[C@@H](C(=O)O)N",
        # "beta-hydroxytyrosine": r"C1=CC(=CC=C1[C@H]([C@@H](C(=O)O)N)O)O",
        # "glutamic acid": r"C(CC(=O)O)[C@@H](C(=O)O)N",
        # "glycine": r"NCC(=O)O",
        # "3,5-dihydroxyphenylglycine": r"N[C@H](C(=O)O)c1cc(O)cc(O)c1",
        # "pipecolic acid": r"C1CCN[C@@H](C1)C(=O)O",
        # "threonine": r"C[C@H]([C@@H](C(=O)O)N)O",
        # "N5-formyl-N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CN(C=O)O",
        # "N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CNO",
        # "ornithine": r"C(C[C@@H](C(=O)O)N)CN",
        # "beta-alanine": r"NCCC(=O)O",
        # "histidine": r"C1=C(NC=N1)C[C@@H](C(=O)O)N",
        # "D-alanine": r"C[C@H](C(=O)O)N",
        # "2-aminoadipic acid": r"C(C[C@@H](C(=O)O)N)CC(=O)O",
        # "lysine": r"C(CCN)C[C@@H](C(=O)O)N",
        # "2-aminoisobutyric acid": r"O=C(O)C(N)(C)C",
        # "2,3-dihydroxybenzoic acid": r"C1=CC(=C(C(=C1)O)O)C(=O)O",
        # "salicylic acid": r"C1=CC=C(C(=C1)C(=O)O)O",
        # "anthranillic acid": r"C1=CC=C(C(=C1)C(=O)O)N",
        # "A", "B", "C", "D",
        ("A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11"): "#039e73",
        ("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11"): "#56b4e9",
        ("C1", "C2", "C4"): "#cc79a7",
        ("D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12"): "#f0e442",
        ("other",): "#808285",
        ("-",): "#ffffff",
    }
    motif_color_map = {name: color for names, color in colors.items() for name in names}

    # Get the length of sequences and number of sequences
    num_sequences = len(sequences)
    sequence_length = len(sequences[0])

    # Create a grid of colors based on the sequences
    color_grid = [
        [motif_color_map.get(str(motif), '#ceccca') for i, motif in enumerate(seq._motifs)]
        for seq in sequences
    ]

    # Convert the color grid into a numpy array for easier handling
    color_grid = np.array(color_grid)

    # Create the plot
    fig, ax = plt.subplots(figsize=(sequence_length + 3, num_sequences + 2))  # Adjust size for names and numbers

    # Create a grid of colored blocks
    for i in range(num_sequences):
        for j in range(sequence_length):
            ax.add_patch(plt.Rectangle((j, num_sequences - i - 1), 1, 1, color=color_grid[i, j], ec='black'))

    # Add names to the left of each row
    for i, name in enumerate([seq._identifier for seq in sequences]):
        ax.text(-0.5, num_sequences - i - 0.5, name, ha='right', va='center', fontsize=20, fontweight='bold')

    # Add column numbers at the bottom
    for j in range(sequence_length):
        ax.text(j + 0.5, -0.5, str(j + 1), ha='center', va='center', fontsize=20, fontweight='bold')

    # Add text inside each cell
    for i in range(num_sequences):
        for j in range(sequence_length):
            # truncate text if it is too long
            if len(str(sequences[i][j])) > 10:
                text = str(sequences[i][j])[:10] + "..."
            else:
                text = str(sequences[i][j])
            if text == "-":
                text = ""
            ax.text(j + 0.5, num_sequences - i - 0.5, text, ha='center', va='center', fontsize=8, color='black', fontweight='bold')

    # Set axis limits
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(0, num_sequences)

    # Turn off the axis
    ax.axis('off')

    # Save the plot as a PNG file
    plt.savefig(output, dpi=300, bbox_inches='tight')
    plt.close(fig)

    pass

def main() -> None:
    args = cli()
    
    # parse out data
    npaids = []
    assignments = []
    primary_sequences_str = []
    with open(args.cluster_assignments, "r") as f:
        f.readline()
        for line in f:
            line = line.strip().split("\t")
            npaid, _, _, primary_sequence, _, _, _, assignment, *_ = line
            npaids.append(npaid)
            assignments.append(assignment)
            primary_sequences_str.append(primary_sequence)
    print(len(npaids), len(assignments), len(primary_sequences_str))

    # create sequence objects from primary sequences
    primary_sequences = []
    for i, sequence_str in enumerate(primary_sequences_str):
        npaid = npaids[i]
        motifs_str = sequence_str.split("|")[1:-1]
        seq = [SequenceMotif(motif_str) for motif_str in motifs_str]
        primary_sequences.append(Sequence(npaid, seq))
    # print(len(primary_sequences))

    # loop over clusters and extract npaids and sequences
    unique_cluster_assignments = set(assignments)
    unique_cluster_assignments = [int(assignment) for assignment in unique_cluster_assignments]
    unique_cluster_assignments = sorted(unique_cluster_assignments)

    # print(max(unique_cluster_assignments))

    # parse scoring matrix
    scoring_matrix = parse_scoring_matrix(args.scoring_matrix)

    # define score func
    def score_func(a, b):
        if isinstance(a, Gap) or isinstance(b, Gap):
            return 4
        elif isinstance(a, Gap) or isinstance(b, Gap):
            return 0
        return scoring_matrix[(str(a.name), str(b.name))]

    # loop over clusters and extract npaids and sequences, and do analysis
    for cluster_assignment in tqdm(unique_cluster_assignments):
        # print(cluster_assignment)
        cluster_indices = [i for i, assignment in enumerate(assignments) if assignment == str(cluster_assignment)]
        cluster_npaids = [npaids[i] for i in cluster_indices]
        cluster_sequences = [primary_sequences[i] for i in cluster_indices]
        # print(cluster_assignment, len(cluster_npaids), len(cluster_sequences))

        try:
            # align sequences
            result = multiple_sequence_alignment(
                seqs=cluster_sequences,
                gap_penalty=5,
                end_gap_penalty=5,
                score_func=score_func,
            )
        except: 
            print(f"Error in cluster {cluster_assignment}")
            continue

        # get max sequence length
        max_length = max([len(seq) for seq in result])
        for seq in result:
            # print(seq)
            if len(seq) < max_length:
                seq._motifs.extend([Gap() for _ in range(max_length - len(seq))])

        # remove columns in alignment that are all gaps
        gap_columns = []
        for i in range(max_length):
            column_items = [seq._motifs[i] for seq in result]
            if all([isinstance(item, Gap) for item in column_items]):
                gap_columns.append(i)
        # print(gap_columns)

        new_result = []
        for seq in result:
            new_seq = [seq._motifs[i] for i in range(max_length) if i not in gap_columns]
            new_result.append(Sequence(seq._identifier, new_seq))
        result = new_result

        # visualize alignment
        out_path = f"{args.output}/alignment_cluster_{cluster_assignment}.png"
        visualize_alignment(result, out_path)

if __name__ == "__main__":
    main()