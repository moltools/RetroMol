#!/usr/bin/env python3
"""
TODO:
    - parse NPAtlas; turn into fingerprints
    - get true compound-BGC pairs; turn into fingerprints
    - plot distribution true pair scores vs. random scores (random scores determined by compounds Tc0.0-0.7 to real and Tc0.7-0.9
    - update naive substitution matrix: custom scores for polyketide units, Tc scores for anything else
    - plot retrieval of true compound (in NPAtlas compounds) for every BGC (split on BGC types all/only PK/only NRP)
    - plot shared biosynthetic UMAP? Have to do some feature selection
"""

import argparse
from pathlib import Path

import numpy as np
from tqdm import tqdm

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
)
from retromol_fingerprint.scoring import cosine_counted


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=str, required=True)
    parser.add_argument("--npatlas-retromol", required=True)
    parser.add_argument("--scoring-matrix", required=True)
    return parser.parse_args()


def main() -> None:
    args = cli()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    submat = parse_substitution_matrix(args.scoring_matrix, unknown_token="unknown")
    mod_embs = build_module_embeddings(submat, embedding_dim=4)

    # Test target for retrieval
    target_smi = r"CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"
    target_submission = Submission(name="erythromycin", smiles=target_smi)
    ruleset = RuleSet.load_default(match_stereochemistry=True)  # NOTE: make sure NPAtlas is parsed with same default ruleset! otherwise retrieval won't work properly
    target_result = run_retromol(submission=target_submission, rules=ruleset)
    target_graph = retromol_result_to_graph(target_result, node_name_attr="name", unidentified_node_name="unknown")
    target_fp = module_graph_fingerprint(graph=target_graph, embeddings=mod_embs, n_bits=2048, node_name_attr="name", counted=True)

    # closest_smi = r"CCC1OC(=O)C(C)C(OC2CC(C)(OC)C(O)C(C)O2)C(C)C(OC2OC(C)CC(N(C)C)C2O)C(C)(O)CC(C)C(=O)C(C)C(O)C1C"
    # closest_submission = Submission(name="closest_to_erythromycin", smiles=closest_smi)
    # ruleset = RuleSet.load_default(match_stereochemistry=False)
    # closest_result = run_retromol(submission=closest_submission, rules=ruleset)
    # closest_graph = retromol_result_to_graph(closest_result, node_name_attr="name", unidentified_node_name="unknown")
    # closest_fp = module_graph_fingerprint(graph=closest_graph, embeddings=mod_embs, n_bits=2048, node_name_attr="name", counted=True)
    # print(cosine_counted(target_fp, target_fp))
    # print(cosine_counted(target_fp, closest_fp))
    # exit("CHK")

    # Parse NPAtlas
    npatlas_fps, npatlas_smis = [], []
    for i, serialized_result in tqdm(enumerate(iter_json(args.npatlas_retromol, jsonl=True))):
        result = Result.from_dict(serialized_result)
        npatlas_smi = result.submission.smiles
        npatlas_graph = retromol_result_to_graph(result, node_name_attr="name", unidentified_node_name="unknown")
        npatlas_fp = module_graph_fingerprint(graph=npatlas_graph, embeddings=mod_embs, n_bits=2048, node_name_attr="name", counted=True)
        npatlas_fps.append(npatlas_fp)
        npatlas_smis.append(npatlas_smi)

        # TODO: delete after testing!
        # if i == 999:
        #     break

    fps = np.stack(npatlas_fps, axis=0)
    print(fps.shape)

    # Compute similarity scores
    scores = np.array([
        cosine_counted(fp, target_fp)
        for fp in npatlas_fps
    ])

    # Top 10 indices (highest scores first)
    top_k = 10
    top_indices = np.argsort(scores)[::-1][:top_k]
    print("\nTop hits:")
    for rank, idx in enumerate(top_indices, start=1):
        print(
            f"{rank:2d}. "
            f"score={scores[idx]:.4f} "
            f"smiles={npatlas_smis[idx]}"
        )



if __name__ == "__main__":
    main()
