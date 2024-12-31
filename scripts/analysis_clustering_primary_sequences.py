#!/usr/bin/env python3
import argparse 
import os
import hashlib
import pickle

import matplotlib.pyplot as plt
from retromol.data.motifs import MOTIFS
from retromol.chem import mol_to_fingerprint
from rdkit import Chem, DataStructs, RDLogger
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import numpy as np
from tqdm import tqdm

from versalign.motif import Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence

from analysis_primary_sequences import parse_results_folder

from sklearn.cluster import KMeans, MiniBatchKMeans
from sklearn.cluster import AgglomerativeClustering
from collections import defaultdict

from scipy.cluster.hierarchy import dendrogram, linkage



# use Arial font
plt.rcParams["font.family"] = "Arial"


motif_dict = {}
for motif in MOTIFS:
    smiles = motif["smiles"]
    mol = Chem.MolFromSmiles(smiles)
    fp = mol_to_fingerprint(mol)
    fp_arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, fp_arr)
    motif_dict[motif["name"]] = fp_arr


def parse_dataset(path: str):
    id_to_classes = {}
    id_to_fp = {}

    with open(path, "r") as f:
        f.readline()
        for line in f:
            line = line.strip().split(",")
            accession, smiles, classes, *_ = line
            try:
                mol = Chem.MolFromSmiles(smiles)
                for atom in mol.GetAtoms():
                    atom.SetIsotope(0)
                fp = mol_to_fingerprint(mol, 2, 2048)
                fp_arr = np.zeros((1,))
                DataStructs.ConvertToNumpyArray(fp, fp_arr)
            except:
                continue
            id_to_fp[accession] = fp_arr
            id_to_classes[accession] = classes

    return id_to_classes, id_to_fp



def encode_fingerprints_numpy(fingerprints):
    # Ensure all inputs are valid binary NumPy arrays
    if not all(isinstance(fp, np.ndarray) and np.all(np.isin(fp, [0, 1])) for fp in fingerprints):
        raise ValueError("All fingerprints must be NumPy arrays containing only 0 and 1.")
    
    # Ensure all fingerprints are of the same size
    fingerprint_size = len(fingerprints[0])
    if not all(len(fp) == fingerprint_size for fp in fingerprints):
        raise ValueError("All fingerprints must have the same length.")
    
    # Ensure all fingerprints are integers
    fingerprints = [fp.astype(int) for fp in fingerprints]
    
    # Initialize the aggregated fingerprint with an array of zeros of the same size
    aggregated_fingerprint = np.zeros(fingerprint_size, dtype=int)

    for i, binary_fingerprint in enumerate(fingerprints):
        # Encode position: Convert position index to binary array of the same size
        position_binary = np.array(
            list(map(int, bin(i)[2:].zfill(fingerprint_size))), dtype=int
        )
        # XOR the fingerprint with the position encoding
        position_encoded = np.bitwise_xor(binary_fingerprint, position_binary)
        # Aggregate using XOR
        aggregated_fingerprint = np.bitwise_xor(aggregated_fingerprint, position_encoded)

    return aggregated_fingerprint


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="path to retromol results folder")
    parser.add_argument("-a", required=True, type=str, help="path to original retromol input file with annotations")
    parser.add_argument("-o", required=True, type=str, help="path to output folder")
    parser.add_argument("-c", required=False, type=str, help="path to cluster centroids")
    parser.add_argument("-m", required=False, type=str, help="path to primary sequence fingerprints")
    return parser.parse_args()


def calc_tanimoto_similarity(arr1, arr2):
    """
    Calculates the Tanimoto similarity for binary vectors.
    """
    # Dot product gives the intersection (A ∩ B)
    intersection = np.dot(arr1, arr2)
    # Sum gives the total "1s" in each vector
    sum_arr1 = np.sum(arr1)
    sum_arr2 = np.sum(arr2)
    # Calculate Tanimoto similarity
    tanimoto_sim = intersection / (sum_arr1 + sum_arr2 - intersection)
    return tanimoto_sim


def motif_code_to_fp(centroids, motif_code) -> np.ndarray:
    fps = []
    for motif in motif_code:
        if motif["identity"] in motif_dict:
            fps.append(motif_dict[motif["identity"]])
        else:
            mol = Chem.MolFromSmiles(motif["smiles"])
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
            fp = mol_to_fingerprint(mol)
            fp_arr = np.zeros((1,))
            DataStructs.ConvertToNumpyArray(fp, fp_arr)
            motif_dict[motif["identity"]] = fp_arr
            fps.append(fp_arr)
    
    # loop over fps and get 2-mers out
    kmers = []
    for a, b in zip(fps[:-1], fps[1:]):
        encoded = encode_fingerprints_numpy([a, b])
        kmers.append(encoded)

    # get size of centroids as that is the size of the fingerprint
    fp = np.zeros((centroids.shape[0],))
    # for every 2mer in dictyostatin, find the closest centroid, and add 1 to that bin
    # get 1st nearest neighbor with centroid
    for signature in kmers:
        distances = np.linalg.norm(centroids - signature, axis=1)
        closest = np.argmin(distances)
        fp[closest] = 1
    
    return fp


def parse_motif_code_sequences(centroids, sequences, id_to_class):
    labels = []
    classes = []
    fps = []
    
    for accession, motif_codes in tqdm(sequences.items()):
        if len(motif_codes): # skip compounds without motif codes
            for node, motif_codes in motif_codes.items():
                for motif_code in motif_codes.values():
                    fp = motif_code_to_fp(centroids, motif_code)
                    labels.append(accession)
                    classes.append(id_to_class[accession])
                    fps.append(fp)

    labels = np.array(labels)
    classes = np.array(classes)
    fps = np.array(fps)

    return labels, classes, fps


# 1. Orange: 0xe69f00
# 2. Light blue: 0x56b4e9
# 3. Green: 0x039e73
# 4. Yellow: 0xf0e442
# 5. Dark blue: 0x0072b2
# 6. Red-ish: 0xd55f00
# 7. Pink: 0xcc79a7

# class to color
class_to_color = {
    "Antimycins": "#e69f00"
}


def plot_biosynthetic_space(ax, fps, classes, labels, dim_red_type="TSNE", annotated=True):

    # reformat classes, get out split and sort as list and then make string again
    classes = [c.split("|") for c in classes]
    classes = ["|".join(sorted(c)) for c in classes]
    classes = np.array(classes)

    # PCA
    if dim_red_type == "PCA":
        pca = PCA(n_components=2)
        pcs = pca.fit_transform(fps)
        ev = pca.explained_variance_ratio_
        print(f"PC1 ({ev[0]:.2f})")
        print(f"PC2 ({ev[1]:.2f})")
    # t-SNE
    elif dim_red_type == "TSNE":
        distance_metric = "cosine"
        tsne = TSNE(n_components=2, random_state=0, perplexity=100, n_iter=2000, metric=distance_metric)
        pcs = tsne.fit_transform(fps)
    else:
        raise ValueError("Invalid dimensionality reduction type.")

    # plot
    if annotated is True:
        for i, c in enumerate(np.unique(classes)):
            idx = classes == c
            # use rainbow color palette
            color = class_to_color.get(c, plt.cm.rainbow(i / len(np.unique(classes))))
            ax.scatter(pcs[idx, 0], pcs[idx, 1], label=c, alpha=0.5, s=10, edgecolors="k", linewidth=0.5, color=color)
        ax.legend()

        # ALSO DRAW annotateed plot with plotly so that i can interactively see the labels
        import plotly.express as px
        import pandas as pd

        df = pd.DataFrame(pcs, columns=["PC1", "PC2"])
        df["class"] = classes
        df["label"] = labels

        fig = px.scatter(df, x="PC1", y="PC2", color="class", hover_data=["label"])
        fig.show()

    else:
        # draw all dot sin gray
        ax.scatter(pcs[:, 0], pcs[:, 1], alpha=0.5, s=10, edgecolors="k", linewidth=0.5, color="#ceccca")


    # remove ticks and tick labels
    ax.tick_params(axis="both", which="both", bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)

    # remove lines around plot, the box
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["left"].set_visible(False)


def cluster_fingerprints(fps, classes, og_labels, output_path):
    # only keep items with where class contains "Open-chain polyketides"
    idx = np.array(["Lipopeptides" in c for c in classes])
    fps = fps[idx]
    classes = classes[idx]
    og_labels = og_labels[idx]
    print(fps.shape, classes.shape, og_labels.shape)

    # only keep fingerprints where at least 3 bits are set
    idx = np.sum(fps, axis=1) >= 3
    fps = fps[idx]
    classes = classes[idx]
    print(fps.shape, classes.shape, og_labels.shape)

    # cluster fingerprints using hierarchical clustering
    n_clusters = 50
    clustering = AgglomerativeClustering(n_clusters=n_clusters, linkage="ward")
    labels = clustering.fit_predict(fps)

    cluster_indices = defaultdict(list)
    for idx, cluster_id in enumerate(labels):
        cluster_indices[cluster_id].append(og_labels[idx])
    
    # Print indices per cluster
    for cluster_id, cluster_items in cluster_indices.items():
        print(f"Cluster {cluster_id}: {cluster_items}")

def main():
    args = cli()

    RDLogger.DisableLog("rdApp.*")

    # parse arguments
    results_path = args.i
    dataset_path = args.a
    output_path = args.o
    centroids_file = args.c
    motif_info = args.m

    if not centroids_file:
        # get all 2-mer fingerprints
        kmers = set()
        for a, motif_a in enumerate(MOTIFS):
            for b, motif_b in enumerate(MOTIFS):
                kmers.add((motif_a["name"], motif_b["name"]))
        kmers = list(kmers)
        kmer_fps = []
        for (a, b) in tqdm(kmers):
            fp_a = motif_dict[a]
            fp_b = motif_dict[b]
            fps = [fp_a, fp_b]
            encoded = encode_fingerprints_numpy(fps)
            kmer_fps.append(encoded)
        kmer_fps = np.array(kmer_fps)

        kmeans = MiniBatchKMeans(n_clusters=2048, random_state=0, n_init="auto").fit(kmer_fps)
        centroids = kmeans.cluster_centers_
        binary_centroids = (centroids > 0.5).astype(int)
        # save binary centroids to numpy file
        np.save(os.path.join(output_path, "cluster_centroids.npy"), binary_centroids)
    else:
        binary_centroids = np.load(centroids_file)
    print("loaded centroids:", binary_centroids.shape)

    # example similarity metric 1
    dictyostatin = Chem.MolFromSmiles(r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O")
    discodermolide = Chem.MolFromSmiles(r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O")
    tanimito_sim = DataStructs.TanimotoSimilarity(mol_to_fingerprint(dictyostatin), mol_to_fingerprint(discodermolide))
    print(tanimito_sim)

    motif_code = ["ethanoic acid", "C1", "C2", "B2", "B2", "D2", "C2", "B2", "C1", "B1", "B2", "B2"]
    motif_code = [{"identity": motif} for motif in motif_code]
    fp1 = motif_code_to_fp(binary_centroids, motif_code)
    motif_code = ["ethanoic acid", "C1", "C2", "B2", "B1", "D2", "D2", "B2", "C1", "B1", "B2", "C1", "C1"]
    motif_code = [{"identity": motif} for motif in motif_code]
    fp2 = motif_code_to_fp(binary_centroids, motif_code)
    biosyn_sim = calc_tanimoto_similarity(fp1, fp2)
    print(biosyn_sim)

    # example similarity metric 2
    megalomycin = Chem.MolFromSmiles(r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O")
    deoxyerythronolide = Chem.MolFromSmiles(r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C")
    tanimito_sim = DataStructs.TanimotoSimilarity(mol_to_fingerprint(megalomycin), mol_to_fingerprint(deoxyerythronolide))
    print(tanimito_sim)

    motif_code = ["ethanoic acid", "B6", "B2", "A2", "D6", "B2", "B2"]
    motif_code = [{"identity": motif} for motif in motif_code]
    fp1 = motif_code_to_fp(binary_centroids, motif_code)
    motif_code = ["ethanoic acid", "B2", "B2", "A2", "D2", "B2", "B2"]
    motif_code = [{"identity": motif} for motif in motif_code]
    fp2 = motif_code_to_fp(binary_centroids, motif_code)
    biosyn_sim = calc_tanimoto_similarity(fp1, fp2)
    print(biosyn_sim)

    # parse inputs
    id_to_class, id_to_fp = parse_dataset(dataset_path)
    _, _, _, sequences = parse_results_folder(results_path)
    print(len(id_to_class), len(id_to_fp), len(sequences))

    # parse motif code sequences
    if not motif_info:
        labels, classes, fps = parse_motif_code_sequences(binary_centroids, sequences, id_to_class)
        data_tuple = (labels, classes, fps)
        # save all three objects to single numpy file, use pickle
        with open(os.path.join(output_path, "fingerprints_motif_codes.pkl"), "wb") as f:
            pickle.dump(data_tuple, f)
    else:
        with open(motif_info, "rb") as f:
            labels, classes, fps = pickle.load(f)
    print(fps.shape, classes.shape, labels.shape)

    # create figure
    fig = plt.figure(figsize=(8, 8))
    gs = fig.add_gridspec(2, 4)

    # create subplots
    ax1 = fig.add_subplot(gs[0:2, 0:4])  # tanimoto similarity
    # plot_biosynthetic_space(ax1, fps, classes, labels, annotated=True)

    # cluster biosynthetic fingerprints
    cluster_fingerprints(fps, classes, labels, output_path)

    # write out plot
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "clustering_primary_sequences.png"), dpi=300)
    plt.savefig(os.path.join(output_path, "clustering_primary_sequences.svg"))

    

if __name__ == "__main__":
    main()
