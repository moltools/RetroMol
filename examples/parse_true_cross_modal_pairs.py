#!/usr/bin/env python3

import argparse
import pickle
import os
import random
from pathlib import Path

import yaml
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from tqdm import tqdm
from rdkit import RDLogger
from umap import UMAP

from retromol.model.rules import RuleSet, ReactionRule, MatchingRule
from retromol.model.result import Result
from retromol.io.streaming import run_retromol_stream, stream_table_rows
from retromol_fingerprint.fingerprint import module_graph_fingerprint
from retromol_antismash.utils import download_and_prepare
from retromol_antismash.io import AntiSmashOptions, load_regions
from retromol_antismash.inference.model_paras import ParasModel
from retromol_antismash.inference.registry import register_domain_model, annotate_region
from retromol_antismash.modules import NRPSModule, linear_readout


RDLogger.DisableLog('rdApp.*')


NUM_BITS = 2048


ORANGE_HEX = r"#e69f00"
GRAY_HEX = r"#ceccca"
RED_HEX = r"#d55f00"
GREEN_HEX = r"#039e73"
BLUE_HEX = r"#56b4e9"
YELLOW_HEX = r"#f0e442"
PINK_HEX = r"#cc79a7"


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, required=True, help="Path to true cross-modal pairs.")
    parser.add_argument("--cache", type=str, required=True, help="Path to cache dir for downloads.")
    parser.add_argument("--paras", type=str, required=True, help="Path to PARAS model.")
    parser.add_argument("--rxn", type=str, required=True, help="Path to reaction rules.")
    parser.add_argument("--mxn", type=str, required=True, help="Path to matching rules.")
    parser.add_argument("--scoring", type=str, required=True, help="Path to scoring matrix.")
    parser.add_argument("--output", type=str, required=True, help="Path to output PNG.")
    parser.add_argument("--mibig", type=str, required=False, help="Path to dir with MIBiG GBKs")
    return parser.parse_args()


def read_scoring_matrix(path: Path | str) -> pd.DataFrame:
    """
    Reads a scoring matrix from a file. The file can be in either Excel or CSV format.

    :param path: The path to the scoring matrix file.
    :return: A pandas DataFrame containing the scoring matrix.
    """
    # Could be xlsx or csv
    if path.endswith(".xlsx"):
        return pd.read_excel(path)
    elif path.endswith(".csv"):
        return pd.read_csv(path)
    else:
        raise ValueError(f"Unsupported file format: {path}")
    

def load_ruleset(path_rxn: Path | str, path_mxn: Path | str) -> RuleSet:
    """
    Load RetroMol ruleset.

    :param path_rxn: Path to reaction rules YAML file.
    :param path_mxn: Path to matching rules YAML file.
    :return: A RuleSet object containing the loaded rules.
    """
    path_reaction_rules = Path(path_rxn)
    path_matching_rules = Path(path_mxn)

    with open(path_reaction_rules, "r") as fo:
        reaction_rules_data = yaml.safe_load(fo)

    with open(path_matching_rules, "r") as fo:
        matching_rules_data = yaml.safe_load(fo)

    reaction_rules = [ReactionRule.from_dict(d) for d in reaction_rules_data]
    matching_rules = [MatchingRule.from_dict(d) for d in matching_rules_data]

    return RuleSet(False, reaction_rules, matching_rules)


def score_fingerprints(fp1: np.ndarray, fp2: np.ndarray, mode: str) -> float:
    """
    Available modes: cosine
    """
    if mode == "cosine":
        # cosine similarity
        dot_product = np.dot(fp1, fp2)
        norm_fp1 = np.linalg.norm(fp1)
        norm_fp2 = np.linalg.norm(fp2)
        if norm_fp1 == 0 or norm_fp2 == 0:
            return 0.0
        else:
            return dot_product / (norm_fp1 * norm_fp2)
    else:
        raise ValueError(f"Unknown mode: {mode}")


def main() -> None:
    args = cli()

    os.makedirs(args.cache, exist_ok=True)

    scoring_matrix = read_scoring_matrix(args.scoring)

    # randomize column and row headers (in same order) so that every column and row in matrix gets random other column/row name
    scoring_matrix = read_scoring_matrix(args.scoring)

    # labels = [c for c in scoring_matrix.columns if c != "identity"]
    # rng = np.random.default_rng(42)
    # shuffled = rng.permutation(labels)
    # label_map = dict(zip(labels, shuffled))
    # scoring_matrix_shuffled = scoring_matrix.copy()
    # # rename the identity values in the first column
    # scoring_matrix_shuffled["identity"] = scoring_matrix_shuffled["identity"].map(label_map)
    # # rename the score columns
    # scoring_matrix_shuffled = scoring_matrix_shuffled.rename(columns=label_map)
    # scoring_matrix = scoring_matrix_shuffled

    # loop over rows of scoring matrix and convert to dict of dicts
    # columnd and row names are the same and correspond to monomer identities
    substitution_matrix = {}
    for _, row in scoring_matrix.iterrows():
        iden = row["identity"]
        substitution_matrix[iden] = {other_iden: row[other_iden] for other_iden in scoring_matrix.columns if other_iden != "identity"}

    # Add UNIDENTIFIED identity with score 0 to all other identities
    substitution_matrix["UNIDENTIFIED"] = {other_iden: 0 for other_iden in substitution_matrix.keys()}
    for iden in substitution_matrix.keys():
        substitution_matrix[iden]["UNIDENTIFIED"] = 0
    substitution_matrix["UNIDENTIFIED"]["UNIDENTIFIED"] = 0


    # path parsed compounds pickle
    path_parsed_compounds = os.path.join(args.cache, "parsed_compounds.pkl")
    if not os.path.exists(path_parsed_compounds):

        # Parse compounds
        source_iter = stream_table_rows(args.data, sep=",", chunksize=20_000)
        ruleset = load_ruleset(args.rxn, args.mxn)
        print(ruleset)

        pbar = tqdm()

        fprs_compounds = {}
        for evt in run_retromol_stream(
            ruleset=ruleset,
            row_iter=source_iter,
            smiles_col="smiles",
            workers=7,
            batch_size=2_000,
            pool_chunksize=50,
            maxtasksperchild=2_000,
        ):
            # We ignore any errors
            if evt.error is not None:
                pass
            elif evt.result is not None:
                # Result is serialized
                result = Result.from_dict(evt.result)
                cov = result.calculate_coverage()
                if cov <= 0.5:
                    continue  # skip low coverage results
                
                assembly_graph = result.linear_readout.assembly_graph
                G = assembly_graph.g
                for n in G.nodes:
                    molnode = G.nodes[n]["molnode"]
                    if molnode is not None:
                        iden = molnode.identity.name if molnode.identified else "UNIDENTIFIED"
                    else:
                        iden = "UNIDENTIFIED"
                    # Assign iden to node_name_attr
                    G.nodes[n]["name"] = iden

                fpr = module_graph_fingerprint(
                    graph=G,
                    substitution_matrix=substitution_matrix,
                    n_bits=NUM_BITS,
                    radius=2,
                    node_name_attr="name",
                    embedding_dim=16,
                    count_similar_edges=True,
                    counted=True,
                )
                fprs_compounds[result.submission.props["mibig"]] = fpr
                
            else:
                pass
            
            # Progress
            pbar.update(1)

        # Check lengths
        print(len(fprs_compounds))

        # pickle
        with open(path_parsed_compounds, "wb") as fo:
            pickle.dump(fprs_compounds, fo)
    else:
        with open(path_parsed_compounds, "rb") as fo:
            fprs_compounds = pickle.load(fo)
    
    # Download/retrieve relevant GBKs
    if args.mibig is None:
        mibig4_url = r"https://dl.secondarymetabolites.org/mibig/mibig_gbk_4.0.tar.gz"
        mibig4_dir = download_and_prepare(mibig4_url, force=True, cache_dir=args.cache)
    else:
        mibig4_dir = args.mibig

    # make paras prediction to name mapping from matching rules
    mapping = {"unknown": "UNIDENTIFIED"}
    with open(args.mxn, "r") as fo:
        mxn = yaml.safe_load(fo)
    for x in mxn:
        new_name = x["name"]
        collapsed = x["props"]["collapsed"]
        for c in collapsed:
            mapping[c] = new_name

    # loop greedily over files in mibig4_dir and only keep baths where any label is in path name
    gbk_paths = {}
    for path in Path(mibig4_dir).rglob("*.gbk"):
        for label in set(fprs_compounds.keys()):
            if label is not None and label in path.name:
                gbk_paths[label] = path
    print(f"Found {len(gbk_paths)} GBK files corresponding to labels in the data.")

    # Parse GBKs
    path_parsed_gbks = os.path.join(args.cache, "parsed_gbks.pkl")
    if not os.path.exists(path_parsed_gbks):
        pm = ParasModel(threshold=0.1, keep_top=3, cache_dir=args.cache, model_path=args.paras)
        register_domain_model(pm)
        options = AntiSmashOptions(readout_level="cand_cluster")  # or cand_cluster
        fprs_gbks = {}
        for gbk_path in tqdm(gbk_paths.values()):
            file_name = Path(gbk_path).name.split(".")[0]
            compound_fpr = fprs_compounds[file_name]  # should be in there
            regions = load_regions(gbk_path, options)

            best_score = float("-inf")  # find best scoring candidate cluster
            best_fpr = None  # find best scoring candidate cluster
            for region in regions:
                annotate_region(region)
                readout = linear_readout(region)
                modules = readout.biosynthetic_order(by_orf=False)
                modules = [m for m in modules if isinstance(m, NRPSModule)]
                module_names = []
                for module in modules:
                    pred_name = module.predicted_substrate.name
                    new_name = mapping[pred_name]
                    assert new_name in substitution_matrix, f"{new_name} not in substitution_matrix!"
                    module_names.append(new_name)
                
                # make linear graph from module names; store in node name as "name"
                G = nx.Graph()
                for i, name in enumerate(module_names):
                    G.add_node(i, name=name)
                for i in range(len(module_names) - 1):
                    G.add_edge(i, i+1)
                fpr = module_graph_fingerprint(
                    graph=G,
                    substitution_matrix=substitution_matrix,
                    n_bits=NUM_BITS,
                    radius=2,
                    node_name_attr="name",
                    embedding_dim=16,
                    count_similar_edges=True,
                    counted=True,
                )
                score = score_fingerprints(compound_fpr, fpr, mode="cosine")
                if score > best_score:
                    best_score = score
                    best_fpr = fpr

            if best_fpr is not None:
                fprs_gbks[file_name] = best_fpr    

        # pickle
        with open(path_parsed_gbks, "wb") as fo:
            pickle.dump(fprs_gbks, fo)
    else:
        with open(path_parsed_gbks, "rb") as fo:
            fprs_gbks = pickle.load(fo)        

    # Calculate distribution scores true pairs and random pairs
    # get union of keys from fprs_compounds and fprs_gbks and get two matrices of fingerprints, one for compounds and one for gbks, where rows correspond to keys in the union and are ordered the same in both matrices; if a key is missing in one of the dicts it is excluded from the matrices
    keys = set(fprs_compounds.keys()) & set(fprs_gbks.keys())
    print(f"Found {len(keys)} pairs of compounds and GBKs with fingerprints.")
    compound_matrix = np.array([fprs_compounds[k] for k in keys])
    gbk_matrix = np.array([fprs_gbks[k] for k in keys])
    print(f"Compound matrix shape: {compound_matrix.shape}")
    print(f"GBK matrix shape: {gbk_matrix.shape}")

    # # filter out any column in either matrix where column in gbks is all zeros but not in compounds
    # nonzero_gbk_cols = np.where(gbk_matrix.sum(axis=0) > 0)[0]
    # compound_matrix = compound_matrix[:, nonzero_gbk_cols]
    # gbk_matrix = gbk_matrix[:, nonzero_gbk_cols]
    # print(f"After filtering zero columns, compound matrix shape: {compound_matrix.shape}")
    # print(f"After filtering zero columns, GBK matrix shape: {gbk_matrix.shape}")

    # Calculate cosine similarity scores for true pairs
    true_scores = []
    for i in range(len(keys)):
        score = score_fingerprints(compound_matrix[i], gbk_matrix[i], mode="cosine")
        true_scores.append(score)

    # Calcualte cosines similarity scores for random non-true pairs
    num_samples = 1000
    random_scores = []
    for _ in tqdm(range(num_samples)):
        i = np.random.randint(0, len(keys))
        j = np.random.randint(0, len(keys))
        if i != j:
            score = score_fingerprints(compound_matrix[i], gbk_matrix[j], mode="cosine")
            random_scores.append(score)

    # Visualize: figure with 2 panels on left density function of true and random scores, on right PCA on both fingerprints visualized in scatter plot
    plt.rcParams["font.family"] = "Helvetica" 
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), gridspec_kw={"wspace": 0.25})

    # Density plot of scores as smooth filled lines
    bins = np.linspace(0, 1, 21)  # 20 bins from 0 to 1

    true_hist, true_edges = np.histogram(true_scores, bins=bins, density=True)
    random_hist, random_edges = np.histogram(random_scores, bins=bins, density=True)

    true_centers = 0.5 * (true_edges[:-1] + true_edges[1:])
    random_centers = 0.5 * (random_edges[:-1] + random_edges[1:])

    # make smoother x-grid
    x_smooth = np.linspace(0, 1, 300)

    # simple linear interpolation
    true_smooth = np.interp(x_smooth, true_centers, true_hist)
    random_smooth = np.interp(x_smooth, random_centers, random_hist)

    axes[0].plot(
        x_smooth,
        random_smooth,
        label="random pairs",
        color=GRAY_HEX,
        linewidth=3,
    )
    axes[0].fill_between(
        x_smooth,
        random_smooth,
        0,
        color=GRAY_HEX,
        alpha=0.4,
    )

    axes[0].plot(
        x_smooth,
        true_smooth,
        label=f"true pairs (n={len(true_scores)})",
        color=PINK_HEX,
        linewidth=3,
    )
    axes[0].fill_between(
        x_smooth,
        true_smooth,
        0,
        color=PINK_HEX,
        alpha=0.25,
    )

    axes[0].set_xlabel("similarity (cosine)")
    axes[0].set_ylabel("density")
    axes[0].set_yticks([])
    axes[0].set_title("cross-modal similarity scoring")
    axes[0].set_xlim(0, 1)
    axes[0].legend()

    # UMAP scatter plot
    combined_matrix = np.vstack([compound_matrix, gbk_matrix])

    umap = UMAP(
        n_components=2,
        n_neighbors=15,
        min_dist=0.1,
        metric="cosine",
        random_state=42,
    )

    umap_result = umap.fit_transform(combined_matrix)

    compound_umap = umap_result[:len(keys)]
    gbk_umap = umap_result[len(keys):]

    # give little black outline
    axes[1].scatter(compound_umap[:, 0], compound_umap[:, 1], label="products", alpha=0.7, edgecolor="k", linewidth=1, color=GREEN_HEX)
    axes[1].scatter(gbk_umap[:, 0], gbk_umap[:, 1], label="gene clusters", alpha=0.7, edgecolor="k", linewidth=1, color=BLUE_HEX)

    axes[1].set_xlabel("UMAP 1")
    axes[1].set_ylabel("UMAP 2")
    axes[1].set_title("UMAP of biosynthetic fingerprints")

    # no ticks or numbers
    axes[1].set_xticks([])
    axes[1].set_yticks([])
    axes[1].legend()

    plt.savefig(args.output, dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":
    main()
