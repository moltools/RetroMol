import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse 
from populate_database import connect_to_database
from retrieval import SequenceMotif
from versalign.motif import Motif, Gap
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence
from versalign.msa import multiple_sequence_alignment
from collections import Counter
from tqdm import tqdm

plt.rcParams['font.sans-serif'] = "Arial"

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=str, required=True, help="path to output svg file")
parser.add_argument("--scoring-matrix", type=str, required=True, help="path to scoring matrix csv file")
args = parser.parse_args()
out_path = args.output
scoring_matrix_path = args.scoring_matrix

def parse_scoring_matrix(path_to_file: str) -> dict:
    scoring_matrix = {}
    df = pd.read_csv(path_to_file, index_col=0)
    for name in df.columns:
        for other_name, score in df[name].items():
            scoring_matrix[(name, other_name)] = score
    return scoring_matrix
scoring_matrix = parse_scoring_matrix(scoring_matrix_path)

cur, conn = connect_to_database()
query_name = "BGC0000336"
query_sequence = ['tryptophan', 'asparagine', 'aspartic acid', 'threonine', 'glycine', 'ornithine', 'aspartic acid', 'alanine', 'aspartic acid', 'glycine', 'serine', 'glutamic acid', 'phenylalanine']
compound_ids = [
    'NPA020427', 'NPA003812', 'NPA000264', 'NPA020422', 'NPA020420', 
    'NPA007316', 'NPA020421', 'NPA015353', 'NPA018911', 'NPA020425', 
    'NPA001560', 'NPA013266', 'NPA000662', 'NPA019734', 'NPA011480', 
    'NPA020270', 'NPA020426', 'NPA001700', 'NPA020424', 'NPA020423', 
    'NPA020317'
]

query = f"""
    SELECT 
        c.compound_id, 
        cps.primary_sequence_id, 
        STRING_AGG(m.motif_name, '|' ORDER BY psm.position) AS primary_sequence 
    FROM compounds c 
    JOIN compounds_primary_sequences cps ON c.compound_id = cps.compound_id 
    JOIN primary_sequences_motifs psm ON cps.primary_sequence_id = psm.primary_sequence_id 
    JOIN motifs m ON psm.motif_id = m.motif_name 
    WHERE c.compound_id IN %s
    GROUP BY c.compound_id, cps.primary_sequence_id;
"""
cur.execute(query, (tuple(compound_ids),))

def visualize_alignment(sequences: list, output: str) -> None:
    # Define a mapping of residues to colors
    colors = {
        ("A", "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11"): "#039e73",
        ("B", "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11"): "#56b4e9",
        ("C", "C1", "C2", "C4"): "#cc79a7",
        ("D", "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12"): "#f0e442",
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
            # ax.add_patch(plt.Rectangle((j, num_sequences - i - 1), 1, 1, color=color_grid[i, j], ec='black'))
            ax.add_patch(plt.Rectangle((j, num_sequences - i - 1), 1, 1, color=color_grid[i, j]))

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
    plt.savefig(output)
    plt.close(fig)

# define score func
def score_func(a, b):
    if isinstance(a, Gap) or isinstance(b, Gap): return 4
    elif isinstance(a, Gap) or isinstance(b, Gap): return -8
    return scoring_matrix[(str(a.name), str(b.name))]

seqs = set()
for compound_id, _, primary_sequence in tqdm(cur.fetchall()):
    seqs.add((compound_id, primary_sequence))
print(f"fetched {len(seqs)} sequences")

new_seqs = []
for compound_id, primary_sequence in seqs:
    seq = Sequence(compound_id, [SequenceMotif(m) for m in primary_sequence.split("|")[1:-1]])
    new_seqs.append(seq)
print(f"fetched {len(new_seqs)} sequences")
seqs = new_seqs

new_seqs.append(Sequence(query_name, [SequenceMotif(m) for m in query_sequence]))
print(f"added query sequence: {query_name}")
print(f"there are now {len(new_seqs)} sequences")

# align sequences
seqs = multiple_sequence_alignment(
    seqs=seqs,
    gap_penalty=5,
    end_gap_penalty=7,
    score_func=score_func,
)

# get max sequence length
max_length = max([len(seq) for seq in seqs])
for seq in seqs:
    if len(seq) < max_length:
        seq._motifs.extend([Gap() for _ in range(max_length - len(seq))])

# remove columns in alignment that are all gaps
gap_columns = []
for i in range(max_length):
    column_items = [seq._motifs[i] for seq in seqs]
    if all([isinstance(item, Gap) for item in column_items]):
        gap_columns.append(i)

new_seqs = []
for seq in seqs:
    new_seq = [seq._motifs[i] for i in range(max_length) if i not in gap_columns]
    new_seqs.append(Sequence(seq._identifier, new_seq))
seqs = new_seqs

print(f"removed {len(gap_columns)} gap columns")

visualize_alignment(seqs, out_path)
