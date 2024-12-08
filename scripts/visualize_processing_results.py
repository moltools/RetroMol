import argparse 
import os
import json

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# use Arial font in plots
plt.rcParams["font.family"] = "Arial"


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, required=True, help="path to results")
    parser.add_argument("-o", "--output", type=str, required=True, help="path to out folder")
    return parser.parse_args()


def plot_coverage_scores(ax, data_folder):
    succeeded = 0
    failed = 0
    coverage_scores = []
    # read csv file
    processing_results_file = os.path.join(data_folder, "processing_results.csv")
    with open(processing_results_file, "r") as f:
        f.readline()  # skip header
        for line in f:
            _, is_success, coverage_score, *_ = line.strip().split(",")
            if len(coverage_score) > 0:
                coverage_scores.append(float(coverage_score))
            if is_success == "succeeded":
                succeeded += 1
            else:
                failed += 1

    print(len(coverage_scores))
    print(f"success: {succeeded}")
    print(f"failed: {failed}")


    # create barplot for coverage scores
    bins = np.arange(0, 102, 1) 
    hist_values, bin_edges = np.histogram(coverage_scores, bins=bins)
    cumulative_counts = np.cumsum(hist_values)

    ax.bar(bin_edges[:-1], hist_values, width=1, color="#56b4e9", edgecolor="black", linewidth=0.7, zorder=2)
    ax.set_xlabel("heavy atoms of input in identifiable motifs", fontsize=12)
    ax.set_ylabel("compounds per bin (#)", fontsize=12)

    percentiles = [50, 75]
    percentile_values = np.percentile(coverage_scores, percentiles)
    split_value1 = percentile_values[0]
    split_value2 = percentile_values[1]

    # normalize cumulative counts to max bin count
    cumulative_counts = cumulative_counts / max(cumulative_counts) * max(hist_values)
    ax.plot(bin_edges[:-1], cumulative_counts, color="k", linestyle="-", linewidth=2, label="cumulative count", zorder=1)
    ax.fill_between(bin_edges[:-1], cumulative_counts, where=(bin_edges[:-1] <= split_value1), color="#d55f00", alpha=0.3, zorder=0)
    ax.fill_between(bin_edges[:-1], cumulative_counts, where=((bin_edges[:-1] > split_value1) & (bin_edges[:-1] <= split_value2)), color="#f0e442", alpha=0.3, zorder=0)
    ax.fill_between(bin_edges[:-1], cumulative_counts, where=(bin_edges[:-1] > split_value2), color="#039e73", alpha=0.3, zorder=0)

    # set x-axis title
    ax.set_xlim(-1, 101)
    ax.set_xticks(np.arange(0, 101, 20))
    ax.set_xticklabels([f"{i}%" for i in range(0, 101, 20)], fontsize=12)
    ax.set_yticks(np.arange(0, max(hist_values)+1, 100))
    ax.set_yticklabels([f"{i}" for i in range(0, max(hist_values)+1, 100)], fontsize=12)

    # little bit of space on x-axis left and right
    ax.set_xlim(-3, 103)

    # set lines for every y-axis tick
    for i in range(0, max(hist_values)+1, 100):
        ax.axhline(i, color="k", linestyle="--", linewidth=0.5, zorder=1, alpha=0.3)

    ax.text(44, 910, f"cumulative (%)", ha="center", va="bottom", fontsize=12, color="k", rotation=0, zorder=3)

    for i, (p, val) in enumerate(zip(percentiles, percentile_values)):
        correction = 0.0
        if i == 0: correction = 0.3
        if i == 1: correction = -0.2
        ax.axvline(val+correction, color="k", linestyle="--", linewidth=2, zorder=3, alpha=0.8)
        perc_str = "$^\mathrm{th}$%"
        ax.text(val+3, max(hist_values)*0.275, f"{p}{perc_str}", ha="left", va="bottom", fontsize=12, color="k", rotation=0, zorder=3)


def plot_number_sequenceable_nodes(ax, data_folder):
    num_backbones = []
    # folder contains folders, list all these folder paths that start with "restuls_"
    results_folders = [os.path.join(data_folder, f) for f in os.listdir(data_folder) if f.startswith("results_")]
    for folder in results_folders:
        # in that folder there is a file called out.json
        out_file = os.path.join(folder, "out.json")
        # check if file exists
        if not os.path.exists(out_file):
            continue
        with open(out_file, "r") as f:
            out_data = json.load(f)
            # count num backbones
            num_backbones.append(sum([1 for _, identity in out_data["encoding_to_identity"].items() if identity == "_backbone"]))
    
    # create barplot for number of sequenceable nodes
    max_num_bins = max(num_backbones) + 1
    bins = np.arange(0, max_num_bins, 1)
    hist_values, bin_edges = np.histogram(num_backbones, bins=bins)
    ax.bar(bin_edges[:-1], hist_values, width=1, color="#56b4e9", edgecolor="black", linewidth=1.2, zorder=2)
    ax.set_xlabel("number of sequenceable backbones", fontsize=12)
    ax.set_ylabel("compounds per bin (#)", fontsize=12)
    ax.set_xticks(np.arange(0, max_num_bins - 1, 1))
    ax.set_yticks(np.arange(0, max(hist_values)+1, 200))
    ax.set_xticklabels([f"{i}" for i in range(0, max_num_bins - 1, 1)], fontsize=12)
    ax.set_yticklabels([f"{i}" for i in range(0, max(hist_values)+1, 200)], fontsize=12)
    # set y-axis lines
    for i in range(0, max(hist_values)+1, 200):
        ax.axhline(i, color="k", linestyle="--", linewidth=0.5, zorder=1, alpha=0.3)
    # set counts on bin 4 and bin 5
    # ax.text(4, 20, f"#{hist_values[4]}", ha="center", va="bottom", fontsize=12, color="k", rotation=0, zorder=3)
    # ax.text(5, 20, f"#{hist_values[5]}", ha="center", va="bottom", fontsize=12, color="k", rotation=0, zorder=3)
    

def plot_number_of_sequences_per_sequenceable_node(ax, data_folder):
    num_backbones = []
    # folder contains folders, list all these folder paths that start with "restuls_"
    results_folders = [os.path.join(data_folder, f) for f in os.listdir(data_folder) if f.startswith("results_")]
    for folder in results_folders:
        # in that folder there is a file called out.json
        out_file = os.path.join(folder, "out.json")
        # check if file exists
        if not os.path.exists(out_file):
            continue
        with open(out_file, "r") as f:
            out_data = json.load(f)
            # count num backbones
            num_backbones.append(sum([len(motif_codes) for encoding, motif_codes in out_data["encoding_to_motif_codes"].items() if encoding in out_data["encoding_to_identity"] and out_data["encoding_to_identity"][encoding] == "_backbone"]))
    
    # create barplot for number of sequenceable nodes
    max_num_bins = max(num_backbones) + 1
    bins = np.arange(0, max_num_bins, 1)
    hist_values, bin_edges = np.histogram(num_backbones, bins=bins)
    ax.bar(bin_edges[:-1], hist_values, width=1, color="#56b4e9", edgecolor="black", linewidth=1.2, zorder=2)
    ax.set_xlabel("number of sequences per backbone", fontsize=12)
    ax.set_ylabel("compounds per bin (#)", fontsize=12)
    ax.set_xticks(np.arange(0, max_num_bins, 2))
    ax.set_yticks(np.arange(0, max(hist_values)+1, 200))
    ax.set_xticklabels([f"{i}" for i in range(0, max_num_bins, 2)], fontsize=12)
    ax.set_yticklabels([f"{i}" for i in range(0, max(hist_values)+1, 200)], fontsize=12)
    # set y-axis lines
    for i in range(0, max(hist_values)+1, 200):
        ax.axhline(i, color="k", linestyle="--", linewidth=0.5, zorder=1, alpha=0.3)
    # set counts on bin 4 and bin 5
    # ax.text(4, 20, f"#{hist_values[4]}", ha="center", va="bottom", fontsize=12, color="k", rotation=0, zorder=3)
    # ax.text(5, 20, f"#{hist_values[5]}", ha="center", va="bottom", fontsize=12, color="k", rotation=0, zorder=3)

def main() -> None:
    args = cli()
    data_folder = args.input

    # create plot with multiple subplots, 2 rows, 3 columns
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 3, figure=fig)
    plot_coverage_scores(fig.add_subplot(gs[0, 0]), data_folder)
    plot_number_sequenceable_nodes(fig.add_subplot(gs[0, 1]), data_folder)
    plot_number_of_sequences_per_sequenceable_node(fig.add_subplot(gs[0, 2]), data_folder)

    # save the figure
    plt.tight_layout()
    plt.savefig(f"{args.output}/coverage_scores.png", dpi=300)

    
if __name__ == "__main__":
    main()