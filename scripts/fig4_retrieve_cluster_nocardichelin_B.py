import argparse
from typing import Generator
from pathlib import Path

import numpy as np
import umap
import matplotlib.pyplot as plt
from tqdm import tqdm

from retromol_database.duckdb import FINGERPRINT_SIZE
from retromol_fingerprint.fingerprint import TOKEN_UNK, Vocabulary, Fingerprinter
from retromol_alignment.scoring import create_tanimoto_scoring_matrix, HARDCODED_PK_SCORING
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import Converter
from retromol_alignment.msa import calculate_msa
from retromol_antismash.io import AntiSmashOptions, load_regions
from retromol_antismash.inference.registry import register_domain_model, annotate_region
from retromol_antismash.inference.model_paras import ParasModel
from retromol_antismash.modules import linear_readout, PKSModule, NRPSModule
from retromol.model.rules import RuleSet


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gbk-dir", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    return parser.parse_args()


def modules_to_primary_sequence(modules) -> list[list[str]]:
    per_monomer_tokens = []
    for module in modules:
        match module:
            case PKSModule():
                extender_unit_type = module.substrate.extender_unit.value
                per_monomer_tokens.append(["PK", extender_unit_type])
            case NRPSModule():
                tokens = []
                if substrate := module.substrate:
                    if substrate.name == "unknown":
                        tokens.append(TOKEN_UNK)
                    else:
                        tokens.append(substrate.name)
                else:
                    tokens.append(TOKEN_UNK)
                per_monomer_tokens.append(tokens)
            case _:
                raise ValueError(f"Unknown module {module}")
    return per_monomer_tokens


def gbk_to_primary_sequences(gbk_path: Path, options: AntiSmashOptions) -> Generator[list[list[str]], None, None]:
    regions = load_regions(gbk_path, options=options)
    for region in regions:
        qualifiers_str = ", ".join([q.lower() for q in region.qualifiers.get("product", [])])
        is_siderophore = "siderophore" in qualifiers_str
        annotate_region(region)
        modules = linear_readout(region).modules
        modules = sorted(modules, key=lambda x: x.start, reverse=False)
        primary_sequence = modules_to_primary_sequence(modules)
        if is_siderophore:
            primary_sequence.append(["siderophore"])
        if len(primary_sequence) >= 2:
            yield primary_sequence


def main() -> None:
    args = cli()

    # List all GBK paths gbk-dir
    gbk_paths = [p for p in args.gbk_dir.glob("*.gbk")]
    print(f"Found {len(gbk_paths)} GBK files")

    # Setup fingerprinter
    rules = RuleSet.load_default()
    vocab_tokens = set()
    for rule in rules.matching_rules:
        rule_name = rule.name
        vocab_tokens.add(rule_name)
        rule_pseudonyms = rule.pseudonyms
        vocab_tokens.update(set(rule_pseudonyms))
    vocab = Vocabulary(vocab_tokens)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)

    # Setup aligner
    records = []
    for rule in rules.matching_rules:
        rule_name = rule.name
        rule_smiles = rule.smiles
        records.append((rule_name, rule_smiles))
    scoring_matrix = create_tanimoto_scoring_matrix(
        records,
        radius=2,
        num_bits=2048,
        stereochemistry=False,
        self_score_tokens=[TOKEN_UNK, "A", "B", "C", "D"],
        self_score=1.0,
        hardcoded_scores=HARDCODED_PK_SCORING,
    )
    aligner = setup_aligner(scoring_matrix, mode="global")
    converter = Converter(
        to_identifier=lambda s: scoring_matrix.alphabet.index(s),
        from_identifier=lambda i: scoring_matrix.alphabet.index(i),
    )

    # Setup GBK parser
    cache_dir = args.out_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    options = AntiSmashOptions(
        readout_level="cand_cluster",
    )
    register_domain_model(ParasModel(threshold=0.1, keep_top=1, cache_dir=cache_dir))

    # Fingerprint GBKs
    gbk_labels = []
    gbk_seqs = []
    gbk_fingerprints = []
    for gbk_path in tqdm(gbk_paths):
        for idx, per_monomer_tokens in enumerate(gbk_to_primary_sequences(gbk_path, options)):
            fingerprint = fingerprinter.encode(per_monomer_tokens=per_monomer_tokens)
            gbk_labels.append(f"{gbk_path.name}_cand{idx + 1}")
            gbk_seqs.append(per_monomer_tokens)
            gbk_fingerprints.append(fingerprint)
    gbk_fingerprints = np.array(gbk_fingerprints)
    print(gbk_fingerprints.shape)

    # Add compound fingerprint
    compound_tokens = [
        ["serine"],
        ["salicylic acid"],
        ["N-(5-aminopentyl)hydroxylamine", "siderophore"],
        ["butanedioic acid", "siderophore"],
        ["N-(5-aminopentyl)hydroxylamine", "siderophore"],
    ]
    compound_fingerprint = [fingerprinter.encode(per_monomer_tokens=compound_tokens)]
    compound_fingerprint = np.array(compound_fingerprint)
    print(compound_fingerprint.shape)
    print(compound_fingerprint.sum())

    fingerprints = np.concatenate([gbk_fingerprints, compound_fingerprint], axis=0)
    print(fingerprints.shape)

    # Visualize 2D PCDA of GBK fingerprint space
    reducer = umap.UMAP(n_components=2, random_state=42, n_neighbors=15, min_dist=0.5, metric="cosine")
    coords = reducer.fit_transform(fingerprints)
    print(coords.shape)

    gbk_coords = coords[:gbk_fingerprints.shape[0], :]
    compound_coords = coords[gbk_fingerprints.shape[0], :]
    print(compound_coords)

    # get index of gbk fingerprint with label "CP108381.region022.gbk_cand2"
    loi = "CP108381.region022.gbk_cand2"
    loi_idx = gbk_labels.index(loi)
    loi_coords = gbk_coords[loi_idx, :]

    plt.figure(figsize=(6, 6))
    plt.scatter(gbk_coords[:, 0], gbk_coords[:, 1], s=50, alpha=0.8, edgecolor="k", linewidths=1, marker="o", color="#56b4e9")
    plt.scatter(compound_coords[0], compound_coords[1], s=50, alpha=1.0, edgecolor="k", linewidths=1, marker="s", color="#039e73")
    plt.scatter(loi_coords[0], loi_coords[1], s=50, alpha=1.0, edgecolor="k", linewidths=1, marker="o", color="#039e73")
    plt.xlabel(f"UMAP 1")
    plt.ylabel(f"UMAP 2")
    plt.xticks([])
    plt.yticks([])
    plt.tight_layout()
    plt.savefig(args.out_dir / "gbk_space.png", dpi=300)
    plt.savefig(args.out_dir / "gbk_space.svg")
    plt.close()

    # List closest GBK neighbors by cosine distance using only NumPy/math
    query = compound_fingerprint.reshape(-1)  # shape: (n_bits,)
    X = gbk_fingerprints.astype(float)  # shape: (n_gbks, n_bits,)
    q = query.astype(float)
    dot_products = X @ q
    X_norms = np.linalg.norm(X, axis=1)
    q_norm = np.linalg.norm(q)
    cosine_similarities = dot_products / (X_norms * q_norm)
    cosine_distances = 1.0 - cosine_similarities
    k = 20
    closest_indices = np.argsort(cosine_distances)[:k]
    for index in closest_indices:
        dist = cosine_distances[index]
        gbk_label = gbk_labels[index]
        gbk_tokens = gbk_seqs[index]
        print(dist, index, gbk_label, gbk_tokens)


if __name__ == "__main__":
    main()