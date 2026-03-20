#!/usr/bin/env python3

import argparse
from pathlib import Path

import yaml
import pandas as pd
from tqdm import tqdm

from retromol.model.rules import RuleSet, ReactionRule, MatchingRule
from retromol.model.submission import Submission
from retromol.model.result import Result
from retromol.io.streaming import run_retromol_stream, stream_table_rows
from retromol_fingerprint.fingerprint import module_graph_fingerprint


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", type=str, required=True, help="Path to compound dataset.")
    parser.add_argument("--rxn", type=str, required=True, help="Path to reaction rules.")
    parser.add_argument("--mxn", type=str, required=True, help="Path to matching rules.")
    parser.add_argument("--scoring", type=str, required=True, help="Path to scoring matrix.")
    parser.add_argument("--output", type=str, required=True, help="Path to output CSV file with parsed compounds.")
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



def main() -> None:
    """
    Main function to execute the script.
    """
    args = cli()

    source_iter = stream_table_rows(args.data, sep=",", chunksize=20_000)

    scoring_matrix = read_scoring_matrix(args.scoring)
    substitution_matrix = {}
    # loop over rows of scoring matrix and convert to dict of dicts
    # columnd and row names are the same and correspond to monomer identities
    for _, row in scoring_matrix.iterrows():
        iden = row["identity"]
        substitution_matrix[iden] = {other_iden: row[other_iden] for other_iden in scoring_matrix.columns if other_iden != "identity"}

    # Add UNIDENTIFIED identity with score 0 to all other identities
    substitution_matrix["UNIDENTIFIED"] = {other_iden: 0 for other_iden in substitution_matrix.keys()}
    for iden in substitution_matrix.keys():
        substitution_matrix[iden]["UNIDENTIFIED"] = 0
    substitution_matrix["UNIDENTIFIED"]["UNIDENTIFIED"] = 0

    ruleset = load_ruleset(args.rxn, args.mxn)
    print(ruleset)

    # Quick check if every matching rule name is in substitution matrix
    for mr in ruleset.matching_rules:
        assert mr.name in substitution_matrix, f"Matching rule {mr.name} not found in substitution matrix"

    num_bits = 2048
    
    pbar = tqdm()

    with open(args.output, "w") as fo:
        # write header
        fo.write("smiles,coverage," + ",".join([f"b{i}" for i in range(num_bits)]))

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
                smiles = result.submission.smiles
                cov = result.calculate_coverage()
                
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
                    n_bits=num_bits,
                    radius=2,
                    node_name_attr="name",
                    embedding_dim=16,
                    count_similar_edges=True,
                    counted=True,
                )

                # Write line 
                fo.write(f"{smiles},{cov}," + ",".join([f"{x}" for x in fpr]))

            else:
                pass
            
            # Progress
            pbar.update(1)

            # Every 100 updates; flush output file
            if pbar.n % 100 == 0:
                fo.flush()


if __name__ == "__main__":
    main()
