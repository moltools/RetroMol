#!/usr/bin/env python3
import argparse 
import pickle

import numpy as np
import matplotlib.pyplot as plt


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, required=True, help="path to pickle")
    parser.add_argument("-o", type=str, required=True, help="path to output folder")
    return parser.parse_args()


def main():
    args = cli()
    print(args.i)
    print(args.o)

    with open(args.i, "rb") as f:
        record_names, record_smiles_strings, similarity_matrix = pickle.load(f)
    print(len(record_names), len(record_smiles_strings), similarity_matrix.shape)

    # do PCA on the similarity matrix, and visualize with matplotlib
    from sklearn.decomposition import PCA

    pca = PCA(n_components=2)
    pca.fit(similarity_matrix)
    ev = pca.explained_variance_ratio_
    print(ev)

    similarity_matrix_pca = pca.transform(similarity_matrix)
    print(similarity_matrix_pca.shape)

    # cluster and print clusters, also visualize with matplotlib, every cluster different color
    from sklearn.cluster import KMeans

    num_clusters = 200
    kmeans = KMeans(n_clusters=num_clusters, random_state=0, n_init=10)
    kmeans.fit(similarity_matrix_pca)
    labels = kmeans.labels_
    
    # for every cluster, get all there pairwise scores and calculate the average
    for i in range(num_clusters):
        cluster_indices = np.where(labels == i)[0]
        cluster_scores = []
        for j in range(len(cluster_indices)):
            for k in range(j+1, len(cluster_indices)):
                cluster_scores.append(similarity_matrix[cluster_indices[j], cluster_indices[k]])
        # only print if average score is >0.7
        if len(cluster_scores) > 0 and np.mean(cluster_scores) > 0.7:
            # and if number of unique record names is > 1
            if len(set([record_names[cluster_indices[j]] for j in range(len(cluster_indices))])) > 1:
                print(f"average score cluster {i}: {np.mean(cluster_scores)} (num items: {len(cluster_scores)})")

    # print all id->cluster assignment to file:
    with open(args.o + "/cluster_assignment.csv", "w") as f:
        f.write("record_name,cluster\n")
        #print in csv format
        for i in range(len(record_names)):
            f.write(f"{record_names[i]},{labels[i]}\n")
    

    # # draw every cluster separately
    # for i in range(num_clusters):
    #     plt.scatter(similarity_matrix_pca[labels == i, 0], similarity_matrix_pca[labels == i, 1], label=f"cluster {i}")
    # plt.legend()
    # plt.show()




if __name__ == "__main__":
    main()