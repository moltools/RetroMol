#!/usr/bin/env
import argparse
import json
import os
import typing as ty
from tqdm import tqdm
from retrieval import parse_antismash_json


def cli() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--mibig-jsons", type=str, required=True, help="path to folder containing MIBiG JSON files")
    parser.add_argument("--output", type=str, required=True, help="path to output file")
    return parser.parse_args()


def parse_mibig_jsons(path_mibig_jsons: str) -> ty.Generator[ty.Dict[str, ty.Any], None, None]:
    """Parse MIBiG JSON files."""
    file_paths = os.listdir(path_mibig_jsons)
    file_paths = sorted(file_paths)
    for file_path in file_paths:
        if file_path.endswith(".json"):
            try:
                with open(os.path.join(path_mibig_jsons, file_path), "r") as file:
                    accession = file_path.split(".")[0]
                    yield accession, json.load(file)
            except json.JSONDecodeError:
                print(f"Error parsing {file_path}")
                continue


def main() -> None:
    """Parse MIBiG JSON files."""
    args = cli()
    out_file = open(args.output, "w")
    out_file.write("mibig_accession,iteration,primary_sequence\n")
    for accession, record_data in tqdm(parse_mibig_jsons(args.mibig_jsons)):
        queries = parse_antismash_json(record_data, predict_specificities=True)
        for iteration, query in enumerate(queries):
            iteration += 1
            primary_sequence = "|".join(query["primary_sequence"])
            out_file.write(f"{accession},{iteration},'{primary_sequence}'\n")
            out_file.flush()
    out_file.close()
    exit(1)


if __name__ == "__main__":
    main()
