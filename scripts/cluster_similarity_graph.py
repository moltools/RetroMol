#!/usr/bin/env python
import argparse 
import pickle
import os
import json
import typing as ty
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from scipy.sparse.csgraph import laplacian
from rdkit import Chem
from sklearn.metrics.pairwise import cosine_similarity
from scipy.sparse.csgraph import laplacian
from retromol.chem import mol_to_fingerprint
from rdkit.Chem import DataStructs
from scipy.cluster.hierarchy import linkage
from scipy.cluster.hierarchy import linkage, cut_tree, to_tree
from cluster_sequences import SequenceMotif

# use arial
plt.rcParams['font.sans-serif'] = "Arial"

def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--similarity-data", type=str, required=True, help="path to similarity data pickle")
    parser.add_argument("--parsed-data", type=str, required=True, help="path to parsed data pickle")
    parser.add_argument("--dataset", type=str, required=True, help="npatlas compounds with classes, original input")
    parser.add_argument("--output", type=str, required=True, help="path to output folder")
    parser.add_argument("--npatlas", type=str, required=True, help="path to npatlas data")
    return parser.parse_args()

def mol_to_fp(mol):
    fp = mol_to_fingerprint(mol)
    fp_arr = np.zeros((1,))
    DataStructs.ConvertToNumpyArray(fp, fp_arr)
    return fp_arr

def calc_tanimoto_similarity(arr1, arr2):
    intersection = np.dot(arr1, arr2)
    sum_arr1 = np.sum(arr1)
    sum_arr2 = np.sum(arr2)
    tanimoto_sim = intersection / (sum_arr1 + sum_arr2 - intersection)
    return tanimoto_sim

def load_similarity_data(similarity_data: str) -> tuple:
    with open(similarity_data, "rb") as f: data = pickle.load(f)
    return data

def to_newick(linkage_matrix: np.ndarray, labels: ty.List[str]) -> str:
    tree = to_tree(linkage_matrix, rd=False)
    def get_newick(node, newick, parentdist, leaf_names):
        if node.is_leaf(): return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
        else:
            if len(newick) > 0: newick = "):%.2f%s" % (parentdist - node.dist, newick)
            else: newick = ");"
            newick = get_newick(node.get_left(), newick, node.dist, leaf_names)
            newick = get_newick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
            newick = "(%s" % (newick)
            return newick
    newick = get_newick(tree, "", tree.dist, labels)
    return newick

def color_strips(leafs, colors, labels):
    src = """DATASET_COLORSTRIP
SEPARATOR SPACE
DATASET_LABEL label1
COLOR #ff0000
COLOR_BRANCHES 0
DATA
#9606 #ff0000 Human
"""
    for leaf, color, label in zip(leafs, colors, labels):
        src += f"{leaf} {color} {label}\n"
    return src

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
    print(f"Number of npaids in dataset: {len(npaid_to_class)}")

    # load npatlas data json
    npaid_to_organism = {}
    with open(args.npatlas, "r") as f:
        npatlas_data = json.load(f)
        for record in npatlas_data:
            npaid = record["npaid"]
            genus = record["origin_organism"]["genus"]
            npaid_to_organism[npaid] = genus
    print(f"Number of npaids in npatlas: {len(npaid_to_organism)}")

    # count most common organism
    from collections import Counter

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
    record_sequences = ["|".join([str(m) for m in s._motifs]) for s in record_sequences]
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

    # drawm dendrogram
    distance = 1 - alignment_scores
    condensed_distance = distance[np.triu_indices(distance.shape[0], k=1)]
    linkage_matrix = linkage(condensed_distance, method="ward")

    # save linkage matrix as newick
    record_names = [f"{i}_{record_name}" for i, record_name in enumerate(record_names)]
    newick = to_newick(linkage_matrix, record_names)
    with open(os.path.join(args.output, "dendrogram_first.nwk"), "w") as f:
        f.write(newick)

    genus_counter = Counter([npaid_to_organism[npaid.split("_")[1]] for npaid in record_names])
    print(genus_counter.most_common(10))

    # for every leaf determine if genus is Streptomyces or not, if yes color is red, else blue
    leaf_colors = []
    labels = []
    for record_name in record_names:
        genus = npaid_to_organism[record_name.split("_")[1]]
        if genus == "Streptomyces":
            leaf_colors.append("#e69f00")
        elif genus == "Bacillus":
            leaf_colors.append("#56b4e9")
        elif genus == "Pseudomonas":
            leaf_colors.append("#039e73")
        elif genus == "Microcystis":
            leaf_colors.append("#f0e442")
        elif genus == "Lyngbya":
            leaf_colors.append("#0072b2")
        elif genus == "Sorangium":
            leaf_colors.append("#d55f00")
        elif genus == "Trichoderma":
            leaf_colors.append("#cc79a7")
        else:
            leaf_colors.append("#ceccca")
        labels.append(genus)

    # save color strips
    color_strip_src = color_strips(record_names, leaf_colors, labels)
    with open(os.path.join(args.output, "color_strip.txt"), "w") as f:
        f.write(color_strip_src)


    # exit()

    # cut tree for clusters
    # num_clusters = 500
    # clusters = cut_tree(linkage_matrix, n_clusters=num_clusters)
    # print(clusters.shape)
    # cut under a certain threshold
    threshold = 0.5
    clusters = cut_tree(linkage_matrix, height=threshold)
    num_clusters = len(np.unique(clusters))
    print(f"Number of clusters: {num_clusters}")
    
    # # for every cluster gather all class labels, and all genus labels, afterwards, check for each cluster if they are enriched in a certain class or genus
    # cluster_to_class = {}
    # cluster_to_genus = {}
    # for i, record_name in enumerate(record_names):
    #     record_name = record_name.split("_")[1] # appended with index to give unique name for iTOL
    #     cluster = clusters[i].item()
    #     if cluster not in cluster_to_class:
    #         cluster_to_class[cluster] = []
    #     if cluster not in cluster_to_genus:
    #         cluster_to_genus[cluster] = []
    #     class_ = npaid_to_class[record_name].split("|")
    #     genus = npaid_to_organism[record_name] if record_name in npaid_to_organism else "Unknown"
    #     cluster_to_class[cluster].extend(class_)
    #     cluster_to_genus[cluster].append(genus)

    # # for every cluster, do enrichment analysis, which is a simple Chi-squared test
    # from scipy.stats import chi2_contingency, fisher_exact

    # cluster_to_class_enrichment = {}
    # cluster_to_genus_enrichment = {}
    # for cluster in tqdm(range(num_clusters), desc="Calculating enrichment"):
    #     class_counts = {}
    #     genus_counts = {}
    #     for class_ in cluster_to_class[cluster]:
    #         if class_ not in class_counts:
    #             class_counts[class_] = 0
    #         class_counts[class_] += 1
    #     for genus in cluster_to_genus[cluster]:
    #         if genus not in genus_counts:
    #             genus_counts[genus] = 0
    #         genus_counts[genus] += 1
    #     # calculate enrichment
    #     class_enrichment = chi2_contingency(class_counts)
    #     genus_enrichment = chi2_contingency(genus_counts)
    #     print(class_enrichment, genus_enrichment)
    #     cluster_to_class_enrichment[cluster] = class_enrichment
    #     cluster_to_genus_enrichment[cluster] = genus_enrichment
    # print("Enrichment analysis done")


    # for every cluster, get the smiles strings and calculate average tanimoto similarity, with std
    cluster_to_tanimoto_similarity = {}
    for cluster in tqdm(range(num_clusters), desc="Calculating average Tanimoto similarity"):
        cluster_smiles = []
        for i, cluster_assignment in enumerate(clusters):
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
        
    with open(os.path.join(args.output, "cluster_assignment.tsv"), "w") as f:
        f.write("npaid\tsmiles\tnum_primary_sequences\tprimary_sequence\tavg_tanimoto_similarity\tstd_avg_tanimoto_similarity\tbiosynthetic_class\tgenus\tcluster\n")
        for i, (record_name, record_smiles_string) in enumerate(zip(record_names, record_smiles_strings)):
            record_name = record_name.split("_")[1] # appended with index to give unique name for iTOL
            mol = Chem.MolFromSmiles(record_smiles_string)
            for atom in mol.GetAtoms():
                atom.SetIsotope(0)
            smiles = Chem.MolToSmiles(mol)
            primary_sequence = record_sequences[i]
            num_primary_sequences = record_name_counts[record_name]
            top_assigned_cluster = clusters[i].item()
            avg_similarity, std_similarity = cluster_to_tanimoto_similarity[top_assigned_cluster]
            genus = npaid_to_organism[record_name] if record_name in npaid_to_organism else "Unknown"
            # print(genus)
            f.write(f"{record_name}\t{smiles}\t{num_primary_sequences}\t{primary_sequence}\t{avg_similarity}\t{std_similarity}\t{npaid_to_class[record_name]}\t{genus}\t{top_assigned_cluster}\n")
    print("Cluster assignment saved to cluster_assignment.tsv")



if __name__ == "__main__":
    main()
