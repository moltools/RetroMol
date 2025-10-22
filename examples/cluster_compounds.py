# -*- coding: utf-8 -*-

"""Example script that shows how to cluster retrobiosynthetic fingerprints of multiple compounds using RetroMol."""

from typing import List, Tuple

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

from retromol.api import run_retromol_with_timeout
from retromol.fingerprint import (
    FingerprintGenerator,
    NameSimilarityConfig,
    cosine_similarity,
    polyketide_family_of
)
from retromol.helpers import iter_json
from retromol.io import Input as RetroMolInput, Result
from retromol.rules import get_path_default_matching_rules


COMPOUNDS = [
    ("nocardichelin_B", r"CCCCCCCCCCC/C=C\C(=O)N(CCCCCNC(=O)CCC(=O)N(CCCCCNC(=O)[C@@H]1COC(=N1)C2=CC=CC=C2O)O)O"),
    ("desferrioxamin", r"CC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCNC(=O)CCC(=O)N(O)CCCCCN"),
    ("erythromycin_C", r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"),
    ("megalomycin_A", r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)O)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O[C@H]4C[C@H]([C@H]([C@@H](O4)C)O)N(C)C)C)C)O)(C)O"),
    ("6-deoxyerytrhonolide", r"CC[C@@H]1[C@@H]([C@@H]([C@H](C(=O)[C@@H](C[C@@H]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O)C)O)C)C)C)O)C"),
    ("daptomycin", r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"),
    ("discodermolide", r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O"),
    ("dictyostatin", r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O"),
    ("anthracimycin", r"C[C@@H]1/C=C\C=C\[C@H](OC(=O)[C@@H](C(=O)/C=C(/[C@H]2[C@@H]1C=C[C@@H]3[C@@H]2CC=C(C3)C)\O)C)C"),
    ("chlorotonil", r"C[C@@H]1/C=C\C=C\[C@@H](OC(=O)[C@H](C(=O)C(C(=O)[C@@H]2[C@H]1C=C[C@H]3[C@H]2[C@@H](C=C(C3)C)C)(Cl)Cl)C)C"),
    ("avilamycin", r"COC[C@H]1O[C@H]([C@H]([C@H]([C@@H]1O[C@@H]2O[C@@H]([C@@H]([C@@H]([C@H]2O)O[C@H]3C[C@]4(OC5(O[C@@H]4[C@H](O3)C)C[C@H]([C@@H]([C@H](O5)C)O[C@H]6C[C@H]([C@@H]([C@H](O6)C)OC(c7c(OC)c(Cl)c(O)c(Cl)c7C)=O)O)O)C)OC)C)O)OC)O[C@@H]8OC[C@@H]9O[C@@]%10(O[C@@H]([C@@](C(O)C)([C@@H]%11OCO[C@H]%11%10)O)C)O[C@H]9[C@H]8OC(C(C)C)=O")
]


def similarity_matrix(fps: NDArray[np.float32]) -> NDArray[np.float32]:
    """
    Compute pairwise similarity between binary fingerprints.

    :param fps: 2D numpy array of shape (n_samples, n_features) with binary fingerprints.
    :return: 2D numpy array of shape (n_samples, n_samples) with pairwise cosine similarities.
    """
    n = fps.shape[0]
    sim = np.zeros((n, n), dtype=float)
    for i in range(n):
        a = fps[i]
        for j in range(i, n):
            b = fps[j]
            sim_ij = cosine_similarity(a, b)
            sim[i, j] = sim[j, i] = sim_ij
    return sim


def main() -> None:
    """Main function to align two compounds."""

    # Setup fingerprint generator
    path_default_matching_rules = get_path_default_matching_rules()
    collapse_by_name = ["glycosylation", "methylation"]
    cfg = NameSimilarityConfig(family_of=polyketide_family_of, symmetric=True, family_repeat_scale=1)
    generator = FingerprintGenerator(
        matching_rules_yaml=path_default_matching_rules,
        collapse_by_name=collapse_by_name,
        name_similarity=cfg
    )

    # Parse compounds with RetroMol and get fingerprint readouts
    labels, fps = [], []
    for name, smiles in COMPOUNDS:
        input_data = RetroMolInput(cid=name, repr=smiles)
        result = run_retromol_with_timeout(input_data)
        cov = result.best_total_coverage()
        print(f"Coverage for {name}: {cov:.1%}")
        fingerprint = generator.fingerprint_from_result(result, num_bits=512, counted=True)
        # One compound can have multiple readouts / mappings
        num_mappings = fingerprint.shape[0]
        for idx, _ in enumerate(range(num_mappings)):
            labels.append(f"{name}_{idx}")
        fps.append(fingerprint)
    fps_stack = np.vstack(fps)
    print(f"Fingerprint array shape: {fps_stack.shape}")

    # Compute similarity matrix and plot dendrogram
    cosine_sim_matrix = similarity_matrix(fps_stack)
    print(f"Similarity matrix shape: {cosine_sim_matrix.shape}")
    distance_matrix = 1.0 - cosine_sim_matrix 
    distance_matrix = np.maximum(distance_matrix, 0.0)  # kill tiny negatives
    condensed_distance = squareform(distance_matrix, checks=False)
    Z = linkage(condensed_distance, method="average")
    plt.figure(figsize=(12, 6))
    dendrogram(Z, labels=labels, leaf_rotation=90, leaf_font_size=10)
    plt.ylabel("Cosine distance", fontweight="bold")
    plt.subplots_adjust(bottom=0.35)
    plt.show()


if __name__ == "__main__":
    main()
