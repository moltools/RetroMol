#!/usr/bin/env python
import argparse 
import pickle
import os
import typing as ty
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from tqdm import tqdm
from sklearn.manifold import TSNE
from sklearn.cluster import SpectralClustering
from scipy.sparse.csgraph import laplacian
from rdkit import Chem
from sklearn.metrics.pairwise import cosine_similarity
from scipy.sparse.csgraph import laplacian
from cluster_sequences import SequenceMotif
from retromol.chem import mol_to_fingerprint
from rdkit.Chem import DataStructs
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.cluster.hierarchy import dendrogram, linkage, cut_tree, to_tree

# use arial
plt.rcParams['font.sans-serif'] = "Arial"

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

    return membership_scores, most_probable_cluster, eigenvalues, eigenvectors

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--similarity-data", type=str, required=True, help="path to similarity data pickle")
    parser.add_argument("--parsed-data", type=str, required=True, help="path to parsed data pickle")
    parser.add_argument("--dataset", type=str, required=True, help="npatlas compounds with classes, original input")
    parser.add_argument("--output", type=str, required=True, help="path to output folder")
    return parser.parse_args()

def mol_to_fp(mol):
    fp = mol_to_fingerprint(mol)
    fp_arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, fp_arr)
    return fp_arr


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


def load_similarity_data(similarity_data: str) -> tuple:
    with open(similarity_data, "rb") as f:
        data = pickle.load(f)
    return data

def to_newick(linkage_matrix: np.ndarray, labels: ty.List[str]) -> str:
    """
    Outputs linkage matrix as Newick file.

    Parameters
    ----------
    linkage_matrix : np.ndarray
        condensed distance matrix
    labels : list of str, optional
        leaf labels

    Returns
    -------
    newick : str
        linkage matrix in newick format tree

    Source:
    https://stackoverflow.com/questions/28222179/save-dendrogram-to-newick-format
    """
    tree = to_tree(linkage_matrix, rd=False)

    def get_newick(node, newick, parentdist, leaf_names):
        if node.is_leaf():
            return "%s:%.2f%s" % (
                leaf_names[node.id],
                parentdist - node.dist,
                newick
            )
        else:
            if len(newick) > 0:
                newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else:
                newick = ");"
            newick = get_newick(
                node.get_left(),
                newick,
                node.dist,
                leaf_names
            )
            newick = get_newick(
                node.get_right(),
                ",%s" % (newick),
                node.dist,
                leaf_names
            )
            newick = "(%s" % (newick)
            return newick
    newick = get_newick(tree, "", tree.dist, labels)
    return newick

def main() -> None:
    args = cli()

    # parse original dataset input into npatlas -> class dict
    npaid_to_class = {}
    with open(args.dataset, "r") as f:
        f.readline()
        lines = f.readlines()
        for line in tqdm(lines):
            npaid, smiles, class_ = line.strip().split(",")
            npaid_to_class[npaid] = class_
    # print(npaid_to_class)
    print(f"Number of npaids in dataset: {len(npaid_to_class)}")

    # load pickled data
    record_names, record_smiles_strings, combined_matrix = load_similarity_data(args.similarity_data)
    print(len(record_names), len(record_smiles_strings), combined_matrix.shape)

    # load parsed data
    with open(args.parsed_data, "rb") as f:
        record_names_parsed_data, record_smiles_strings_parsed_data, record_sequences = pickle.load(f)

    # filter out sequences that are too short (<6 units)
    ids_to_keep = []
    for i, sequence in enumerate(record_sequences):
        if len(sequence) >= 6:
            ids_to_keep.append(i)

    record_names_parsed_data = [record_names_parsed_data[i] for i in ids_to_keep]
    record_smiles_strings_parsed_data = [record_smiles_strings_parsed_data[i] for i in ids_to_keep]
    record_sequences = [record_sequences[i] for i in ids_to_keep]
    print(len(record_names_parsed_data), len(record_smiles_strings_parsed_data), len(record_sequences))

    # check if record names are the same
    for i, (record_name, record_name_parsed_data) in enumerate(zip(record_names, record_names_parsed_data)):
        if record_name != record_name_parsed_data:
            print(f"Record names do not match: {record_name} != {record_name_parsed_data}")

    # for every record name, count how many times it appears in the list, create dict with counts
    record_name_counts = {}
    for record_name in record_names:
        if record_name not in record_name_counts:
            record_name_counts[record_name] = 0
        record_name_counts[record_name] += 1

    # reconstruct data
    upper_triangle = np.triu(combined_matrix, k=1)
    lower_triangle = np.tril(combined_matrix, k=-1)
    alignment_scores = upper_triangle + upper_triangle.T
    alignment_lengths = lower_triangle + lower_triangle.T
    alignment_scores = alignment_scores / 100

    # for every pairwise alignment length < 5, set the alignment score to 0
    alignment_scores[alignment_lengths < 5] = 0

    # for every pairwise alignment score < 0.8, set the alignment score to 0
    alignment_scores[alignment_scores < 0.7] = 0

    # if there are multiple sequences with the same sequence, only keep one
    # first gather all sequences per npaid
    npaid_to_sequences = {}
    for i, record_name in enumerate(record_names):
        npaid = record_name
        if npaid not in npaid_to_sequences:
            npaid_to_sequences[npaid] = []
        npaid_to_sequences[npaid].append(i)
    
    # for every npaid, check if there are multiple sequences
    for npaid, sequence_ids in npaid_to_sequences.items():
        if len(sequence_ids) > 1:
            print(sequence_ids)

    exit()

    # per NPAID, dereplicate. If there are multiple sequences, but one has less "other", than pick that one 

    # drawm dendrogram
    distance = 1 - alignment_scores
    condensed_distance = distance[np.triu_indices(distance.shape[0], k=1)]
    linkage_matrix = linkage(condensed_distance, method="ward")

    # save linkage matrix as newick
    record_names = [f"{i}_{record_name}" for i, record_name in enumerate(record_names)]
    newick = to_newick(linkage_matrix, record_names)
    with open(os.path.join(args.output, "dendrogram_first.nwk"), "w") as f:
        f.write(newick)

    # create annotation files
    # 1) length of primary sequences
    # 2) biosynthetic class of primary sequences
    # 3) specief of primary sequences
    # 4) PK/NRP/hybrid of primary sequences


    # plt.figure(figsize=(8, 5))
    # dendrogram(linkage_matrix, labels=record_names, orientation='top', leaf_rotation=90)
    # plt.title('Hierarchical Clustering Dendrogram')
    # plt.xlabel('Sample Index')
    # plt.ylabel('Distance')
    # plt.show()
        
    exit()

    # multiple alignment scores by the alignment lengths
    # alignment_scores = alignment_scores * alignment_lengths

    # peform spectral overlap clustering on the alignment scores
    # cluster, but soft assignment, items can be part of multiple clusters
    num_clusters = 500
    print(f"Performing spectral clustering with {num_clusters} clusters")
    clusters, top_clusters, eigenvalues, eigenvectors = spectral_clustering_with_probable_cluster(alignment_scores, n_clusters=num_clusters)
    # first 100 eigenvalues
    eigenvalues = eigenvalues[:10]
    plt.plot(range(1, len(eigenvalues) + 1), eigenvalues, marker='o')
    plt.xlabel('Index')
    plt.ylabel('Eigenvalue')
    plt.title('Eigenvalue Spectrum')
    plt.savefig(os.path.join(args.output, "eigenvalue_spectrum_first.svg"))
    plt.close()
    print(clusters.shape, top_clusters.shape)

    # visualize clustering in 2D
    # tsne = TSNE(n_components=2, perplexity=50, max_iter=300)
    # print(eigenvectors.shape, eigenvectors[:, 1:3].shape)
    # tsne_embedding = tsne.fit_transform(eigenvectors[:, 1:3])
    # print(tsne_embedding.shape)
    plt.figure(figsize=(10, 10))
    # sns.scatterplot(x=tsne_embedding[:, 0], y=tsne_embedding[:, 1], hue=top_clusters, palette="tab20", legend='full')
    sns.scatterplot(x=eigenvectors[:, 1], y=eigenvectors[:, 2])
    plt.title("Spectral Clustering")
    plt.xlabel("t-SNE 1")
    plt.ylabel("t-SNE 2")
    plt.savefig(os.path.join(args.output, "spectral_clustering_first.svg"))

    # for every cluster, get the smiles strings and calculate average tanimoto similarity, with std
    cluster_to_tanimoto_similarity = {}
    for cluster in tqdm(range(num_clusters), desc="Calculating average Tanimoto similarity"):
        cluster_smiles = []
        for i, cluster_assignment in enumerate(top_clusters):
            if cluster_assignment == cluster:
                cluster_smiles.append(record_smiles_strings[i])
        #format cluster_smiles
        cluster_mols = [Chem.MolFromSmiles(smiles) for smiles in cluster_smiles]
        for mol in cluster_mols:
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
        # get fingerprints
        cluster_fps = [mol_to_fp(mol) for mol in cluster_mols]
        # calculate similarity
        similarities = []
        for i in range(len(cluster_fps)):
            for j in range(i+1, len(cluster_fps)):
                similarity = calc_tanimoto_similarity(cluster_fps[i], cluster_fps[j])
                similarities.append(similarity)
        # calculate average similarity
        avg_similarity = np.mean(similarities)
        std_similarity = np.std(similarities)
        cluster_to_tanimoto_similarity[cluster] = (avg_similarity, std_similarity)
        

    print_top_clusters = 3
    with open(os.path.join(args.output, "cluster_assignment_first.tsv"), "w") as f:
        f.write("npaid\tsmiles\tnum_primary_sequences\tprimary_sequence\tavg_tanimoto_similarity\tstd_avg_tanimoto_similarity\tbiosynthetic_class\t" + '\t'.join(['pick_{0}'.format(i+1) for i in range(print_top_clusters)]) + "\n")
        for i, (record_name, record_smiles_string) in enumerate(zip(record_names, record_smiles_strings)):
            mol = Chem.MolFromSmiles(record_smiles_string)
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
            smiles = Chem.MolToSmiles(mol)
            cluster_assignments = clusters[i]
            top_clusters = np.argsort(cluster_assignments)[::-1][:print_top_clusters]
            primary_sequence = "|".join([str(m) for m in record_sequences[i]._motifs])
            num_primary_sequences = record_name_counts[record_name]
            top_assigned_cluster = top_clusters[0]
            avg_similarity, std_similarity = cluster_to_tanimoto_similarity[top_assigned_cluster]
            f.write(f"{record_name}\t{smiles}\t{num_primary_sequences}\t{primary_sequence}\t{avg_similarity}\t{std_similarity}\t{npaid_to_class[record_name]}\t" + '\t'.join([str(cluster+1) for cluster in top_clusters]) + "\n")
    print("Cluster assignment saved to cluster_assignment_first_clustering.tsv")



if __name__ == "__main__":
    main()
