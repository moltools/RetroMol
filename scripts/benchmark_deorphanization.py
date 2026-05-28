#!/usr/bin/env python3
"""
TODO:
    - plot retrieval of true compound (in NPAtlas compounds) for every BGC (split on BGC types all/only PK/only NRP); with and without alignment-based optimization
    - improve scoring matrix based on alpha/beta amino acid identity, starters, acids etc.
    - move from pyhmmer to original PARAS implementation with hmmer command line tool
    - describe which motifs/monomers get clustered together using singular value decomposition

    - what's up with stereochemsitry substitutation matrix vs. predicted paras substrates???
"""

import argparse
import glob
import os
import json
from pathlib import Path
from collections import Counter, defaultdict

from rdkit import Chem, RDLogger
import networkx as nx
import numpy as np
import pandas as pd
from tqdm import tqdm

import matplotlib.pyplot as plt

from retromol.io.json import iter_json
from retromol.model.result import Result
from retromol.model.submission import Submission
from retromol.model.rules import RuleSet
from retromol.pipelines.parsing import run_retromol
from retromol_fingerprint.fingerprint import (
    retromol_result_to_graph,
    module_graph_fingerprint,
    parse_substitution_matrix,
    build_module_embeddings,
    pks_coarse,
)
from retromol_fingerprint.scoring import cosine_counted
from retromol_antismash.io import AntiSmashOptions, load_regions
from retromol_antismash.inference.registry import DOMAIN_MODELS, GENE_MODELS, register_domain_model, register_gene_model, annotate_region
from retromol_antismash.inference.model_paras import ParasModel
from retromol_antismash.modules import linear_readout, NRPSModule, PKSModule, NRPSSubstrate

from retromol.chem.mol import smiles_to_mol
from retromol.chem.fingerprint import mol_to_morgan_fingerprint, calculate_tanimoto_similarity

from retromol_alignment.scoring import create_scoring_matrix
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import align, Converter

RULESET = RuleSet.load_default(match_stereochemistry=True)


RDLogger.DisableLog("rdApp.*")


def parse_compound(smi: str, mod_embs) -> tuple[np.ndarray, Result]:
    submission = Submission(name="erythromycin", smiles=smi)
    result = run_retromol(submission=submission, rules=RULESET)
    graph = retromol_result_to_graph(result, node_name_attr="name", unidentified_node_name="unknown")
    fp = module_graph_fingerprint(graph=graph, embeddings=mod_embs, node_name_attr="name")
    return fp, result


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=str, required=True)
    parser.add_argument("--npatlas-retromol", required=True)
    parser.add_argument("--scoring-matrix", required=True)
    parser.add_argument("--mibig-gbks", required=True)
    parser.add_argument("--mibig-jsons", required=True)
    return parser.parse_args()


def translate_paras_name(name: str) -> str:
    return {
        "S-adenosylmethionine": "adenosylmethionine",
        "3R-methylglutamic acid": "3-methylglutamic acid",
    }.get(name, name)


def modules_to_graph(modules) -> nx.Graph:
    graph = nx.Graph()

    prev_module_idx = None
    for module_idx, module in enumerate(modules):
        # print(module)
        match module:
            case PKSModule():
                anatomy = module.anatomy
                has_kr = anatomy.has_active_KR
                has_dh = anatomy.has_active_DH
                has_er = anatomy.has_active_ER
                match has_kr, has_dh, has_er:
                    case False, False, False: name = "AX"
                    case True, False, False: name = "BX"
                    case True, True, False: name = "CX"
                    case True, True, True: name = "DX"
                    case _: name = "A"
            case NRPSModule():
                if module.predicted_substrate is not None:
                    predicted: NRPSSubstrate = module.predicted_substrate
                    if predicted is not None:
                        name = translate_paras_name(predicted.name)
                    else:
                        name = "unknown"
                else:
                    name = "unknown"
            case _:
                raise ValueError(f"unknown module type: {type(module)}")

        graph.add_node(module_idx, name=name)
        if prev_module_idx is not None:
            graph.add_edge(prev_module_idx, module_idx)
        prev_module_idx = module_idx

    return graph


def modules_to_list(modules) -> list[str]:
    module_names = []
    for module in modules:
        # print(module)
        match module:
            case PKSModule():
                anatomy = module.anatomy
                has_kr = anatomy.has_active_KR
                has_dh = anatomy.has_active_DH
                has_er = anatomy.has_active_ER
                match has_kr, has_dh, has_er:
                    case False, False, False: name = "AX"
                    case True, False, False: name = "BX"
                    case True, True, False: name = "CX"
                    case True, True, True: name = "DX"
                    case _: name = "A"
            case NRPSModule():
                if module.predicted_substrate is not None:
                    predicted: NRPSSubstrate = module.predicted_substrate
                    if predicted is not None:
                        name = translate_paras_name(predicted.name)
                    else:
                        name = "unknown"
                else:
                    name = "unknown"
            case _:
                raise ValueError(f"unknown module type: {type(module)}")

        module_names.append(name)

    return module_names


score_fps = cosine_counted


def smiles_with_most_heavy_atoms(smiles_list: list[str]) -> str | None:
    best_smiles = None
    best_heavy_atoms = -1

    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue

        num_heavy_atoms = mol.GetNumHeavyAtoms()

        if num_heavy_atoms > best_heavy_atoms:
            best_smiles = smiles
            best_heavy_atoms = num_heavy_atoms

    return best_smiles


def rerank(cluster_modules, top_retromol: list[Result], aligner, converter) -> list[int]:
    cluster_seq = modules_to_list(cluster_modules)

    alignment_scores = []

    for result in top_retromol:
        cumulative_alignment_score = 0
        for path in result.linear_readout.paths:
            compound_seq = []
            for node in path:
                if node.identified:
                    name = pks_coarse(node.identity.name)
                else:
                    name = "unknown"
                compound_seq.append(name)

            score1, _, _ = align(aligner, cluster_seq, compound_seq, converter)
            score2, _, _ = align(aligner, cluster_seq, compound_seq[::-1], converter)
            cumulative_alignment_score += max(score1, score2)

        alignment_scores.append(cumulative_alignment_score)

    # resort indices based on alignment_scores, highest first
    original = list(range(len(top_retromol)))
    reranked = sorted(
        range(len(top_retromol)),
        key=lambda i: alignment_scores[i],
        reverse=True,
    )

    if original != reranked:
        print(f"RERANKED!")

    return reranked


def main() -> None:
    args = cli()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    scoring_matrix_df = pd.read_csv(args.scoring_matrix, sep="\t", index_col=0)
    scoring_matrix_df["unknown"] = 0.0
    scoring_matrix_df.loc["unknown"] = 0.0
    scoring_matrix_df.loc["unknown", "unknown"] = 1.0

    scoring_matrix = create_scoring_matrix(scoring_matrix_df)
    aligner = setup_aligner(
        scoring_matrix,
        mode="global",
        # milder gap penalties
        open_internal_insertion_score=-0.5,
        extend_internal_insertion_score=-0.1,
        open_internal_deletion_score=-0.5,
        extend_internal_deletion_score=-0.1,
        # make terminal gaps less severe/free, extra modules at beginning/end of compoun dpath should not destroy score
        open_left_insertion_score=-0.1,
        extend_left_insertion_score=-0.05,
        open_right_insertion_score=-0.1,
        extend_right_insertion_score=-0.1,
        open_left_deletion_score=-0.1,
        extend_left_deletion_score=-0.05,
        open_right_deletion_score=-0.1,
        extend_right_deletion_score=-0.1,
    )
    converter = Converter(
        to_identifier=lambda s: scoring_matrix.alphabet.index(s),
        from_identifier=lambda i: scoring_matrix.alphabet[i],
    )

    # seq1 = ["BX", "BX", "AX", "DX", "BX", "BX"]
    # seq2 = ["propanoic acid", "BX", "BX", "DX", "BX", "BX"]
    # score, aln1, aln2 = align(aligner, seq1, seq2, converter)
    # print(score)
    # print(aln1)
    # print(aln2)
    #
    # exit("CHK")

    submat = parse_substitution_matrix(args.scoring_matrix, unknown_token="unknown", unknown_score=0.0, self_score=1.0)
    mod_embs = build_module_embeddings(submat, embedding_dim=4)

    # Test target for retrieval
    query1_smi = r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"
    query2_smi = r"CCCCCCCCCC(=O)N[C@@H](CC1=CNC2=CC=CC=C21)C(=O)N[C@H](CC(=O)N)C(=O)N[C@@H](CC(=O)O)C(=O)N[C@H]3[C@H](OC(=O)[C@@H](NC(=O)[C@@H](NC(=O)[C@H](NC(=O)CNC(=O)[C@@H](NC(=O)[C@H](NC(=O)[C@@H](NC(=O)[C@@H](NC(=O)CNC3=O)CCCN)CC(=O)O)C)CC(=O)O)CO)[C@H](C)CC(=O)O)CC(=O)C4=CC=CC=C4N)C"
    query1_fp, _ = parse_compound(smi=query1_smi, mod_embs=mod_embs)
    query2_fp, _ = parse_compound(smi=query2_smi, mod_embs=mod_embs)
    target_fp = None

    options = AntiSmashOptions(readout_level="cand_cluster")

    # first read jsons and get BGC accession + smiles pair
    mibig_pairs = {}
    for mibig_json_path in tqdm(glob.glob(os.path.join(args.mibig_jsons, "*.json"))):
        with open(mibig_json_path, "r") as f:
            data = json.load(f)
        accession = data["accession"]
        # version = data["version"]
        if compounds := data.get("compounds", None):
            structures = []
            for compound in compounds:
                if structure := compound.get("structure", None):
                    structures.append(structure)
            if len(structures):
                structure = smiles_with_most_heavy_atoms(structures)
                mibig_pairs[accession] = structure

    accession_to_gbk_path = {}
    allowed = ["NRPS", "T1PKS", "NRPS;T1PKS"]
    included_all = set()
    included_per_class = defaultdict(list)
    for mibig_gbk_path in tqdm(glob.glob(os.path.join(args.mibig_gbks, "*.gbk"))):
        regions = load_regions(mibig_gbk_path, options=options)
        if not len(regions):
            continue
        lengths = [r.end - r.start for r in regions]
        longest_region_idx = np.argmax(lengths)
        region = regions[longest_region_idx]

        product_classes = region.qualifiers["product"]
        product_classes: list[str] =  sorted(product_classes)
        product_class_str = ";".join(product_classes)
        if product_class_str in allowed and region.id in mibig_pairs:
            accession_to_gbk_path[region.id] = mibig_gbk_path
            included_all.add(region.id)
            included_per_class[product_class_str].append(region.id)
    mibig_pairs = {accession: smiles for accession, smiles in mibig_pairs.items() if accession in included_all}


    # these are all:
    # - non-retired entries
    # - have at least 1 compounda ssociated to it
    # - are exclusively classified as T1PKS, NRPS, or T1PKS+NRPS

    # 660 true pairs:
    # - NRPS+T1PKS: 155
    # - T1PKS: 277
    # - NRPS: 228

    # now we know which mibig accession have at least 1 compound and are in allowed biosynthetic classes
    print(len(mibig_pairs))
    for class_str, accessions in included_per_class.items():
        print(class_str, len(accessions))

    register_domain_model(ParasModel(threshold=0.1, keep_top=1, cache_dir=args.work_dir))

    true_scores_by_class = {}
    random_scores_by_class = {}

    all_fp_pairs = {}
    accession_to_modules = {}

    for allowed_class_str in allowed:
        failed = 0
        succeeded = 0
        true_scores = []
        fp_pairs = []
        pairs = {accession: mibig_pairs[accession] for accession in included_per_class[allowed_class_str]}
        print(allowed_class_str, len(pairs))
        for accession, smiles in tqdm(pairs.items()):
            gbk = accession_to_gbk_path[accession]
            regions = load_regions(gbk, options=options)
            lengths = [r.end - r.start for r in regions]
            longest_region_idx = np.argmax(lengths)
            region = regions[longest_region_idx]
            annotate_region(region)
            try:
                modules = linear_readout(region).modules
                if not len(modules) >= 2:
                    continue
                accession_to_modules[accession] = modules
                graph = modules_to_graph(modules)
                fp_cluster = module_graph_fingerprint(graph=graph, embeddings=mod_embs, node_name_attr="name")

                fp_compound, _ = parse_compound(smi=smiles, mod_embs=mod_embs)

                score = cosine_counted(fp_cluster, fp_compound)
                true_scores.append(score)
                fp_pairs.append((accession, fp_cluster, fp_compound))
                succeeded += 1
            except:
                failed += 1
                continue

        all_fp_pairs[allowed_class_str] = fp_pairs

        random_scores = []
        for _ in range(len(true_scores)):
            # get two random pairs, not the same
            idx1, idx2 = np.random.choice(len(fp_pairs), size=2, replace=False)
            accession1, fp_cluster1, fp_compound1 = fp_pairs[idx1]
            accession2, fp_cluster2, fp_compound2 = fp_pairs[idx2]
            random_scores.append(cosine_counted(fp_cluster1, fp_compound2))

        print(f"failed: {failed} / {succeeded}")
        true_scores_by_class[allowed_class_str] = true_scores
        random_scores_by_class[allowed_class_str] = random_scores

    labels = list(true_scores_by_class.keys())
    true_data = [true_scores_by_class[label] for label in labels]
    random_data = [random_scores_by_class[label] for label in labels]
    x = np.arange(len(labels))
    box_width = 0.25
    pair_spacing = 0.45
    fig, ax = plt.subplots(figsize=(8, 5))
    bp_true = ax.boxplot(
        true_data,
        positions=x - pair_spacing / 2,
        widths=box_width,
        patch_artist=True,
        showmeans=True,
    )
    bp_random = ax.boxplot(
        random_data,
        positions=x + pair_spacing / 2,
        widths=box_width,
        patch_artist=True,
        showmeans=True,
    )
    for box in bp_true["boxes"]:
        box.set_facecolor("tab:blue")
        box.set_alpha(0.5)
    for box in bp_random["boxes"]:
        box.set_facecolor("tab:orange")
        box.set_alpha(0.5)
    ax.set_xticks(x)
    ax.set_xticklabels([l.replace(";", "-") for l in labels])
    # ax.set_xticklabels([l.replace(";", "-") for l in labels], rotation=30, ha="right")
    ax.set_ylim(-0.15, 1.05)
    ax.set_yticks(np.linspace(0, 1, 6))
    ax.set_ylabel("Similarity")
    ax.set_xlabel("BGC class")
    # ax.set_title("True vs random BGC-compound pair scores by class")
    ax.legend(
        [bp_true["boxes"][0], bp_random["boxes"][0]],
        ["True pairs", "Random pairs"],
        loc="lower center",
        bbox_to_anchor=(0.5, 0.02),
        ncol=2,
        frameon=False,
    )
    plt.tight_layout()
    plt.savefig(os.path.join(args.work_dir, "true_vs_random_scores_by_class.png"), dpi=300)
    plt.close()

    # Parse NPAtlas
    npatlas_fps, npatlas_smis, npatlas_retromol = [], [], []
    for i, serialized_result in tqdm(enumerate(iter_json(args.npatlas_retromol, jsonl=True))):
        result = Result.from_dict(serialized_result)
        npatlas_smi = result.submission.smiles
        npatlas_graph = retromol_result_to_graph(result, node_name_attr="name", unidentified_node_name="unknown")
        npatlas_fp = module_graph_fingerprint(graph=npatlas_graph, embeddings=mod_embs, node_name_attr="name")
        npatlas_fps.append(npatlas_fp)
        npatlas_smis.append(npatlas_smi)
        npatlas_retromol.append(result)
        # TODO: delete after testing!
        # if i == 999:
        #     break
    # stack npatlas compounds with true compounds we will try to recover
    for _, smi in tqdm(mibig_pairs.items(), desc="Computing retromol result for mibig compounds"):
        try:
            fp, result = parse_compound(smi=smi, mod_embs=mod_embs)
            npatlas_fps.append(fp)
            npatlas_smis.append(smi)
            npatlas_retromol.append(result)
        except:
            pass
    # fps = np.stack(npatlas_fps, axis=0)

    # loop over all true pairs. for every true pair, we retrieve k closest compounds from npatlas with cluster fp
    retrieval_percentages_by_class = {}
    ks = [1, 5, 10, 50, 100]
    thresholds = [0.4, 0.675]
    for class_str, fp_pairs in tqdm(all_fp_pairs.items()):
        all_top_k_scores = []
        all_top_k_scores_reranked = []

        for accession, cluster_fp, _ in tqdm(fp_pairs):
            associated_smiles = mibig_pairs[accession]
            associated_mol = smiles_to_mol(associated_smiles)
            associated_fp = mol_to_morgan_fingerprint(associated_mol, radius=2, num_bits=2048)
            scores = np.array([score_fps(cluster_fp, fp) for fp in npatlas_fps])
            top_indices = np.argsort(scores)[::-1][:1000]
            top_smis = [npatlas_smis[i] for i in top_indices]
            top_retromol = [npatlas_retromol[i] for i in top_indices]
            cluster_modules = accession_to_modules[accession]
            reranked = rerank(cluster_modules, top_retromol, aligner, converter)
            # reranked_top_indices = [top_indices[i] for i in reranked]

            # default recovery
            top_mols = [smiles_to_mol(s) for s in top_smis]
            top_morgan_fps = [mol_to_morgan_fingerprint(m, radius=2, num_bits=2048) for m in top_mols]

            top_tcs = [calculate_tanimoto_similarity(associated_fp, fp) for fp in top_morgan_fps]
            top_k_scores = []
            for k in ks:
                top_k_scores.append(max(top_tcs[:k]))
            all_top_k_scores.append(top_k_scores)

            # raranked recovery
            reranked_top_tcs = [top_tcs[i] for i in reranked]
            reranked_top_k_scores = []
            for k in ks:
                reranked_top_k_scores.append(max(reranked_top_tcs[:k]))
            all_top_k_scores_reranked.append(reranked_top_k_scores)


        all_top_k_scores = np.array(all_top_k_scores)
        all_top_k_scores_reranked = np.array(all_top_k_scores_reranked)

        retrieval_percentages_by_class[class_str] = {
            "original": {},
            "reranked": {},
        }

        for threshold in thresholds:
            original_percentages = []
            reranked_percentages = []

            for col_idx, k in enumerate(ks):
                original_hits = all_top_k_scores[:, col_idx] >= threshold
                reranked_hits = all_top_k_scores_reranked[:, col_idx] >= threshold

                original_percentages.append(100.0 * original_hits.mean())
                reranked_percentages.append(100.0 * reranked_hits.mean())

            retrieval_percentages_by_class[class_str]["original"][threshold] = original_percentages
            retrieval_percentages_by_class[class_str]["reranked"][threshold] = reranked_percentages

    for class_str, method_to_thresholds in retrieval_percentages_by_class.items():
        fig, ax = plt.subplots(figsize=(6, 4))

        for threshold in thresholds:
            ax.plot(
                ks,
                method_to_thresholds["original"][threshold],
                marker="o",
                linestyle="-",
                label=f"Original, TC ≥ {threshold}",
            )

            ax.plot(
                ks,
                method_to_thresholds["reranked"][threshold],
                marker="o",
                linestyle="--",
                label=f"Reranked, TC ≥ {threshold}",
            )

        ax.set_xscale("log")
        ax.set_xticks(ks)
        ax.set_xticklabels([str(k) for k in ks])

        ax.set_ylim(0, 100)
        ax.set_xlabel("Top-k retrieved compounds")
        ax.set_ylabel("True compound recovered (%)")
        ax.set_title(class_str.replace(";", "-"))
        ax.legend(frameon=False)

        plt.tight_layout()

        safe_class_str = class_str.replace(";", "_").replace("/", "_")
        plt.savefig(
            os.path.join(
                args.work_dir,
                f"retrieval_percentages_{safe_class_str}_reranked.png",
            ),
            dpi=300,
        )
        plt.close()


if __name__ == "__main__":
    main()
