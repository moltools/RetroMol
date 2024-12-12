#!/usr/bin/env python
import argparse 
import pickle
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.cluster import SpectralClustering
from scipy.sparse.csgraph import laplacian
from rdkit import Chem
from sklearn.metrics.pairwise import cosine_similarity
from scipy.sparse.csgraph import laplacian

def spectral_clustering_with_probable_cluster(similarity_matrix, n_clusters=200):
    """
    Perform spectral clustering and compute both soft membership scores and the most probable cluster.

    Parameters:
        similarity_matrix (numpy.ndarray): The similarity matrix (n x n).
        n_clusters (int): Number of clusters to assign.

    Returns:
        tuple:
            - numpy.ndarray: A matrix (n x n_clusters) with membership scores for each item.
            - numpy.ndarray: A vector (n,) indicating the most probable cluster for each item.
    """
    # Step 1: Normalize the similarity matrix to create a graph Laplacian
    laplacian_matrix, diag = laplacian(similarity_matrix, normed=True, return_diag=True)

    # Step 2: Compute the smallest n_clusters + 1 eigenvectors (ignoring the first trivial one)
    eigenvalues, eigenvectors = np.linalg.eigh(laplacian_matrix)
    embedding = eigenvectors[:, 1:n_clusters + 1]  # Use the next n_clusters eigenvectors

    # Step 3: Normalize the embedding (row-wise) for stability
    embedding_normalized = embedding / np.linalg.norm(embedding, axis=1, keepdims=True)

    # Step 4: Compute soft membership scores by cosine similarity to cluster prototypes
    cluster_prototypes = embedding_normalized[:n_clusters]  # Use first n_clusters rows as cluster prototypes
    membership_scores = cosine_similarity(embedding_normalized, cluster_prototypes)

    # Step 5: Normalize membership scores so that they sum to 1 for each item
    membership_scores = membership_scores / membership_scores.sum(axis=1, keepdims=True)

    # Step 6: Determine the most probable cluster for each item
    most_probable_cluster = np.argmax(membership_scores, axis=1)

    return membership_scores, most_probable_cluster

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--similarity-data", type=str, required=True, help="path to similarity data pickle")
    parser.add_argument("--output", type=str, required=True, help="path to output folder")
    return parser.parse_args()

def load_similarity_data(similarity_data: str) -> tuple:
    with open(similarity_data, "rb") as f:
        data = pickle.load(f)
    return data

def matrix_heatmap(matrix: np.ndarray, output: str) -> None:
    plt.figure(figsize=(10, 10))
    sns.heatmap(matrix, cmap="viridis")
    plt.savefig(output)

def main() -> None:
    args = cli()

    # load data
    record_names, record_smiles_strings, combined_matrix = load_similarity_data(args.similarity_data)
    print(len(record_names), len(record_smiles_strings), combined_matrix.shape)

    # reconstruct data
    upper_triangle = np.triu(combined_matrix, k=1)
    lower_triangle = np.tril(combined_matrix, k=-1)
    alignment_scores = upper_triangle + upper_triangle.T
    alignment_lengths = lower_triangle + lower_triangle.T
    alignment_scores = alignment_scores / 100

    # multiple alignment scores by the alignment lengths
    alignment_scores = alignment_scores * alignment_lengths

    # peform spectral overlap clustering on the alignment scores
    # cluster, but soft assignment, items can be part of multiple clusters
    num_clusters = 500
    print(f"Performing spectral clustering with {num_clusters} clusters")
    clusters, top_clusters = spectral_clustering_with_probable_cluster(alignment_scores, n_clusters=num_clusters)
    print(clusters.shape, top_clusters.shape)

    print_top_clusters = 10
    with open(os.path.join(args.output, "cluster_assignment.csv"), "w") as f:
        f.write("record_name,record_smiles_string," + ','.join(['pick_{0}'.format(i + 1) for i in range(print_top_clusters)]) + "\n")
        for i, (record_name, record_smiles_string) in enumerate(zip(record_names, record_smiles_strings)):
            mol = Chem.MolFromSmiles(record_smiles_string)
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
            smiles = Chem.MolToSmiles(mol)
            cluster_assignments = clusters[i]
            top_clusters = np.argsort(cluster_assignments)[::-1][:print_top_clusters]
            f.write(f"{record_name},{smiles}," + ','.join([str(cluster) for cluster in top_clusters]) + "\n")
    print("Cluster assignment saved to cluster_assignment.csv")

    # to t-SNE on the alignment scores and plot
    tsne = TSNE(n_components=2, random_state=0)
    X_2d = tsne.fit_transform(alignment_scores)
    plt.figure(figsize=(6, 5))
    plt.scatter(X_2d[:, 0], X_2d[:, 1], c="#56b4e9", edgecolors="k", s=50, alpha=0.5)
    plt.title("t-SNE on alignment scores")
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.xticks([])
    plt.yticks([])
    plt.savefig(os.path.join(args.output, "t-SNE_alignment_scores.png"), dpi=300)

if __name__ == "__main__":
    main()
