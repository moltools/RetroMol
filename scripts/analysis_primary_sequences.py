#!/usr/bin/env python3
import argparse 
import json
import os
from collections import Counter
from typing import Dict, List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.sankey import Sankey
from rdkit import Chem
from tqdm import tqdm


from retromol.data.motifs import MOTIFS


# use Arial font
plt.rcParams["font.family"] = "Arial"


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="path to retromol results folder")
    parser.add_argument("-a", required=True, type=str, help="path to original retromol input file with annotations")
    parser.add_argument("-o", required=True, type=str, help="path to output folder")
    return parser.parse_args()


def parse_dataset(path: str) -> Dict[str, List[str]]:
    dataset = {}

    with open(path, "r") as f:
        f.readline()
        for line in f:
            line = line.strip().split(",")
            accession, _, classes, *_ = line
            classes = classes.strip('"').split("|")
            dataset[accession] = classes

    return dataset


def parse_results_folder(path: str):
    coverage_total = {}
    coverage_primary_sequences = {}
    successes = {}
    sequences = {}

    # go over subdirs in dir at path and check how many contain an out.json file (is a success)
    for subdir in tqdm(os.listdir(path)):
        subdir_path = os.path.join(path, subdir)
        if not os.path.isdir(subdir_path): continue

        # get accession from subdir name, assign False to successes if no out.json file, else True
        accession = subdir.split("_")[1]
        if not os.path.exists(os.path.join(subdir_path, "out.json")): 
            successes[accession] = False
            continue
        else: successes[accession] = True

        # read out.json file
        with open(os.path.join(subdir_path, "out.json"), "r") as f:
            data = json.load(f)
            
            # get nodes 
            reactant_smiles = data["smiles"]
            encoding_to_smiles = data["encoding_to_smiles"]
            encoding_to_identity = data["encoding_to_identity"]
            encoding_to_motif_codes = data["encoding_to_motif_codes"]

            # get all node ids of nodes that are in encoding_to_identity but not have identity _backbone
            selected_nodes_without_backbones = []  # nodes outside backbones that are identified (as SMILES)
            selected_backbone_nodes = []  # nodes inside backbones that are identified (as SMILES)
            backbone_nodes = []
            for node_id, identity in encoding_to_identity.items():
                if (identity != "_backbone") and (identity is not None): 
                    selected_nodes_without_backbones.append(encoding_to_smiles[node_id])
                elif identity is None:
                    continue
                else: 
                    backbone_nodes.append(node_id)
                    motif_codes = encoding_to_motif_codes[node_id]
                    if len(motif_codes) > 0:
                        if len(motif_codes) == 1: picked_motif_code = motif_codes["0"]
                        else:
                            picked_motif_code = []
                            picked_motif_code_coverage = float("-inf")
                            for _, motif_code in motif_codes.items():
                                smiles_with_identity = [m["smiles"] for m in motif_code if m["identity"] is not None]
                                atom_tags = set()
                                for smiles in smiles_with_identity:
                                    mol = Chem.MolFromSmiles(smiles)
                                    for atom in mol.GetAtoms():
                                        atom_tag = atom.GetIsotope()
                                        if atom_tag != 0: atom_tags.add(atom_tag)
                                if len(atom_tags) > picked_motif_code_coverage:
                                    picked_motif_code = motif_code
                                    picked_motif_code_coverage = len(atom_tags)
                        # get identified nodes from motif_code and add to selected_backbone_nodes
                        for motif in picked_motif_code:
                            if motif["identity"] is not None:
                                selected_backbone_nodes.append(motif["smiles"])

            # get all tags from smiles
            mol = Chem.MolFromSmiles(reactant_smiles)
            all_tags = set()
            for atom in mol.GetAtoms():
                atom_tag = atom.GetIsotope()
                if atom_tag != 0: all_tags.add(atom_tag)

            # get all backbone tags from smiles
            backbone_tags = set()
            for node, identity in encoding_to_identity.items():
                if identity == "_backbone":
                    mol = Chem.MolFromSmiles(encoding_to_smiles[node])
                    for atom in mol.GetAtoms():
                        atom_tag = atom.GetIsotope()
                        if atom_tag != 0: backbone_tags.add(atom_tag)

            # get tags of identified nodes, only backbones
            backbone_tags = set()
            for smiles in selected_backbone_nodes:
                mol = Chem.MolFromSmiles(smiles)
                for atom in mol.GetAtoms():
                    atom_tag = atom.GetIsotope()
                    if atom_tag != 0: backbone_tags.add(atom_tag)

            # get tags of identified nodes, only non-backbones
            non_backbone_tags = set()
            for smiles in selected_nodes_without_backbones:
                mol = Chem.MolFromSmiles(smiles)
                for atom in mol.GetAtoms():
                    atom_tag = atom.GetIsotope()
                    if atom_tag != 0: non_backbone_tags.add(atom_tag)

            # coverage backbone nodes + non-backbone nodes, as percentage of all tags
            coverage = (len(backbone_tags) + len(non_backbone_tags)) / len(all_tags)
            if coverage > 1: print(f"coverage > 1: {coverage}")
            coverage_total[accession] = coverage
            
            # coverage backbone nodes, as percentage of all tags
            try: coverage = len(backbone_tags) / len(backbone_tags)
            except ZeroDivisionError: coverage = 0.0
            if coverage > 1: print(f"coverage > 1: {coverage}")
            coverage_primary_sequences[accession] = coverage

            # add all motif codes for accession to sequences
            sequences[accession] = encoding_to_motif_codes

    
    return coverage_total, coverage_primary_sequences, successes, sequences


def sankey_plot(ax: plt.Axes, dataset: Dict[str, List[str]], successes: Dict[str, bool]) -> None:
    # define constants
    start_num_compounds = 36454

    # create sankey plot
    sankey = Sankey(ax=ax, scale=1e-5)

    # first split: how many compounds included
    compounds_in_dataset = len(dataset)
    compounds_not_in_dataset = start_num_compounds - compounds_in_dataset
    sankey.add(flows=[start_num_compounds, -compounds_not_in_dataset, -compounds_in_dataset], orientations=[0, 0, -1])#, labels=["NPAtlas", "Not in dataset", "In dataset"])

    # second split: how many compouds succeeded
    compounds_success = sum(successes.values())
    compounds_failure = compounds_in_dataset - compounds_success
    sankey.add(flows=[compounds_in_dataset, -compounds_failure, -compounds_success], orientations=[0, 1, 0], prior=0, connect=(2, 0))#, labels=["", "Failed", "Succeeded"])

    # third split: compound classifications (class is class if only one class, other mixed polyketide)
    compound_classes = Counter()
    for _, classes in dataset.items():
        if len(classes) == 1: compound_classes[classes[0]] += 1
        else: 
            new_class_name = " & ".join(classes)
            # compound_classes["Mixed class"] += 1
            compound_classes[new_class_name] += 1

    # get class names and counts, and sort on class name alphabetically
    class_names = list(compound_classes.keys())
    class_counts = list(compound_classes.values())
    class_names, class_counts = zip(*sorted(zip(class_names, class_counts), key=lambda x: x[1], reverse=True))

    # add flow that separates classes
    flows = [compounds_success]
    orientations = [0]
    labels = [""]
    for class_name, class_count in zip(class_names, class_counts):
        flows.append(-class_count)
        orientations.append(1)
        labels.append(class_name)
    sankey.add(flows=flows, orientations=orientations, prior=1, connect=(2, 0), labels=labels)
    
    # remove axis lines
    diagrams = sankey.finish()
    for diagram in diagrams:
        for text in diagram.texts:
            text.set_fontsize(8)
            for text in diagram.texts:
                label = text.get_text()
                if len(label) > 0:
                    text.set_text(label.replace("\n", " "))
                    
    ax.axis("off")


def coverage_plot(ax: plt.Axes, coverage_total: Dict[str, float], coverage_primary_sequences: Dict[str, float]) -> None:
    # create 2d hist for coverages
    coverages_total = list(coverage_total.values())
    coverages_primary_sequences = list(coverage_primary_sequences.values())

    # print values for compounds hat have >90% coverage total but <10% coverage primary sequences
    # for accession, coverage in coverage_total.items():
    #     if coverage > 0.9 and coverage_primary_sequences[accession] < 0.1:
    #         print(accession, coverage, coverage_primary_sequences[accession])

    # coverages to percentages
    coverages_total = [c * 100 for c in coverages_total]
    coverages_primary_sequences = [c * 100 for c in coverages_primary_sequences]

    hist2d = ax.hist2d(coverages_primary_sequences, coverages_total, bins=10, cmap="Reds")
    bin_values = hist2d[0]
   
    # sum values per quadrant
    total = sum(bin_values.flatten())
    top_right = int(sum(bin_values[5:, 5:].flatten()))
    bottom_left = int(sum(bin_values[:5, :5].flatten()))
    top_left = int(sum(bin_values[:5, 5:].flatten()))
    bottom_right = int(sum(bin_values[5:, :5].flatten()))
    print(top_right, bottom_left, top_left, bottom_right)

    ax.set_xlabel("coverage primary sequences")
    ax.set_ylabel("coverage total")
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)

    # make percentages on axes
    ax.set_xticks(range(0, 101, 20))
    ax.set_yticks(range(0, 101, 20))
    ax.set_xticklabels([f"{x}%" for x in range(0, 101, 20)])
    ax.set_yticklabels([f"{y}%" for y in range(0, 101, 20)])

    # put in lines for quadrants and display counts for quadrants
    ax.axvline(50, color="black", linestyle="--", linewidth=0.7)
    ax.axhline(50, color="black", linestyle="--", linewidth=0.7)
    ax.text(30, 70, f"n = {top_left}", fontsize=8)
    ax.text(30, 30, f"n = {bottom_left}", fontsize=8)
    ax.text(60, 70, f"n = {top_right}", fontsize=8)
    ax.text(60, 30, f"n = {bottom_right}", fontsize=8)

    # add color bar
    cbar = plt.colorbar(ax.collections[0], ax=ax)
    # add color bar title
    # cbar.set_label("Number of compounds")


def bar_plot_num_backbones(ax: plt.Axes, sequences):
    # get number of backbones per compound
    num_backbones = []
    print(f"sequences contains motif codes for {len(sequences)} nodes")
    for accession, motif_codes in sequences.items():
        num_backbones.append(len(motif_codes))
        if len(motif_codes) > 5:
            print(accession)

    # sort counts based on number of backbones
    num_backbones_count = Counter(num_backbones)
    num_backbones = sorted(num_backbones)

    # create histogram
    ax.hist(num_backbones, bins=range(0, 11), color="#56b4e9", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("number of backbones")
    ax.set_ylabel("number of compounds")
    max_num_backbones = max(num_backbones)
    max_count = max(Counter(num_backbones).values())
    # make sure x ticks are in the middle of the bars
    tick_positions = [x + 0.5 for x in range(0, max_num_backbones + 1)]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(range(0, max_num_backbones + 1))
    ax.set_xlim(-0.5, max_num_backbones + 1)
    ax.set_ylim(0, max_count + (0.1 * max_count))
    ax.set_yticks(range(0, max_count + 400, 500))

    # put counts in for bars that have height under 100
    for i, count in enumerate(Counter(num_backbones).values()):
        print(i, count)
        if count < 400:
            ax.text(i + 1, count + 10, count, ha="center", va="bottom", fontsize=8)

    # set horizontal lines 
    for height in range(0, max_count + 400, 500):
        ax.axhline(height, color="black", linestyle="--", linewidth=0.7, zorder=0)


def bar_plot_seq_lengths(ax: plt.Axes, sequences):
    # get sequence lengths per compound
    seq_lengths = []
    for accession, motif_codes in sequences.items():
        if len(motif_codes): # skip compounds without motif codes
            for node, motif_codes in motif_codes.items():
                length = max(len(motif_code) for motif_code in motif_codes.values())
                seq_lengths.append(length)

    # sort counts based on sequence length
    seq_lengths = sorted(seq_lengths)
    
    # create histogram
    max_length = max(seq_lengths)
    ax.hist(seq_lengths, bins=range(0, max_length + 1), color="#56b4e9", edgecolor="black", linewidth=0.5)
    ax.set_xlabel("sequence length")
    ax.set_ylabel("number of compounds")
    max_count = max(Counter(seq_lengths).values())
    # make sure x ticks are in the middle of the bars
    tick_positions = [x + 0.5 for x in range(0, max_length + 1) if x % 10 == 0]
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(range(0, max_length + 1, 10))
    ax.set_xlim(-0.5, max_length + 1)
    ax.set_ylim(0, max_count + (0.1 * max_count))
    ax.set_yticks(range(0, max_count + 50, 100))

    # set horizontal lines
    for height in range(0, max_count + 50, 100):
        ax.axhline(height, color="black", linestyle="--", linewidth=0.7, zorder=0)
    

def bar_plot_seq_content(ax: plt.Axes, sequences):
    seq_types = {
        "PK": 0,
        "NRP": 0,
        "hybrid": 0,
        # "other": 0,
    }

    # get counts
    for accession, motif_codes in sequences.items():
        if len(motif_codes): # skip compounds without motif codes
            for node, motif_codes in motif_codes.items():
                for motif_code in motif_codes.values():
                    # check if none of the motifs have a identity None
                    if all(m["identity"] is not None for m in motif_code):
                        identities = set(m["identity"] for m in motif_code)
                        tags = set()
                        for identity in identities:
                            for x in MOTIFS:
                                if x["name"] == identity:
                                    tags.update(x["tags"])
                        if all(tag in tags for tag in ["polyketide", "amino acid"]):
                            seq_types["hybrid"] += 1
                        elif all(tag in tags for tag in ["polyketide"]):
                            seq_types["PK"] += 1
                        elif all(tag in tags for tag in ["amino acid"]):
                            seq_types["NRP"] += 1
                        else:
                            # seq_types["other"] += 1
                            print(f"compound {accession} has other type of sequence")

    # make bar plot with seq_types
    ax.bar(seq_types.keys(), seq_types.values(), color="#56b4e9", edgecolor="black", linewidth=0.5)
    ax.set_ylabel("number of compounds")
    max_count = max(seq_types.values())
    ax.set_ylim(0, max_count + (0.1 * max_count))
    ax.set_yticks(range(0, max_count + 50, 300))
    # rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha="right")

    # set horizontal lines
    for height in range(0, max_count + 50, 300):
        ax.axhline(height, color="black", linestyle="--", linewidth=0.7, zorder=0)


def pmi_plot_seq(ax: plt.Axes, sequences):
    # parse out sequences
    all_sequences = []
    for accession, motif_codes in sequences.items():
        if len(motif_codes): # skip compounds without motif codes
            for node, motif_codes in motif_codes.items():
                for motif_code in motif_codes.values():
                    # check if none of the motifs have a identity None
                    if all(m["identity"] is not None for m in motif_code):
                        sequence = []
                        for motif in motif_code:
                            # find in motif lib
                            motif_tags = []
                            for x in MOTIFS:
                                if x["name"] == motif["identity"]:
                                    motif_tags = x["tags"]
                                    break
                            if "amino acid" in motif_tags:
                                sequence.append("amino acid")
                            elif "polyketide" in motif_tags:
                                sequence.append(motif["identity"])
                            else:
                                sequence.append("other")
                        all_sequences.append(sequence)

    # invert sequences, add start and end token
    all_sequences = [sequence[::-1] for sequence in all_sequences]
    all_sequences = [["start"] + sequence + ["end"] for sequence in all_sequences]

    # print sequences
    # for seq in all_sequences:
    #     print(seq)

    # get bigram counts
    pmi = {}
    for seq in all_sequences:
        for i in range(1, len(seq)):
            bigram = tuple(seq[i-1:i+1])
            if bigram in pmi:
                pmi[bigram] += 1
            else:
                pmi[bigram] = 1
    
    # get unigram counts
    unigram = {}
    for seq in all_sequences:
        for token in seq:
            if token in unigram:
                unigram[token] += 1
            else:
                unigram[token] = 1

    # calculate PMI
    pmi_values = {}
    for bigram, count in pmi.items():
        unigram_1 = bigram[0]
        unigram_2 = bigram[1]
        val = count / (unigram[unigram_1] * unigram[unigram_2])
        # take log2 of val
        pmi_values[bigram] = np.log2(val)

    # sort the labels: <START>, <END>, X, AA, others alphatically
    unigram = dict(sorted(unigram.items(), key=lambda x: x[0]))

    # print all unique tokens
    print(unigram.keys())

    labels = [
        "start",
        "end",
        "other",
        "amino acid",
        "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11",
        "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11",
        "C1", "C2", "C4",
        "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
    ]

    # plot PMI values in heatmap, pos 1 on y_axis, pos 2 on x_axis
    pmi_matrix = []
    for i, token_1 in enumerate(labels):
        row = []
        for j, token_2 in enumerate(labels):
            bigram = (token_2, token_1)
            if bigram in pmi_values:
                row.append(pmi_values[bigram])
            else:
                # add value that is read as NaN by imshow
                row.append(np.nan)
        pmi_matrix.append(row)

    # plot heatmap
    ax.imshow(pmi_matrix, cmap="Reds")
    ax.set_xticks(range(len(labels)))
    ax.set_yticks(range(len(labels)))
    ax.set_xticklabels(labels)
    ax.set_yticklabels(labels, ha="left")

    # rotate x-axis labels
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90, ha="right")
    ax.yaxis.tick_right()

    # add color bar
    # cbar = plt.colorbar(ax.images[0], ax=ax)

    # add grid lines
    ax.set_xticks([x - 0.5 for x in range(1, len(labels))], minor=True)
    ax.set_yticks([y - 0.5 for y in range(1, len(labels))], minor=True)
    ax.grid(which="minor", color="black", linestyle="--", linewidth=0.7)




def main() -> None:
    args = cli()

    # parse arguments
    results_path = args.i
    dataset_path = args.a
    output_path = args.o

    # parse inputs
    dataset = parse_dataset(dataset_path)
    coverage_total, coverage_primary_sequences, successes, sequences = parse_results_folder(results_path)
    print(len(dataset), len(successes), len(coverage_total), len(coverage_primary_sequences), len(sequences))

    # create figure with 3 rows and 4 columns, need to be able to merge cells (use gridspec)
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(4, 4)
    
    # create plots
    ax1 = fig.add_subplot(gs[0:3, 0])
    sankey_plot(ax1, dataset, successes)
    ax2 = fig.add_subplot(gs[3, 0])
    coverage_plot(ax2, coverage_total, coverage_primary_sequences)
    ax3 = fig.add_subplot(gs[3, 1])
    bar_plot_num_backbones(ax3, sequences)
    ax4 = fig.add_subplot(gs[3, 2])
    bar_plot_seq_lengths(ax4, sequences)
    ax5 = fig.add_subplot(gs[3, 3])
    bar_plot_seq_content(ax5, sequences)
    ax6 = fig.add_subplot(gs[0:3, 1:4])
    pmi_plot_seq(ax6, sequences)

    # write out plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "primary_sequence_analysis.png"), dpi=300)
    plt.savefig(os.path.join(output_path, "primary_sequence_analysis.svg"))


if __name__ == "__main__":
    main()
