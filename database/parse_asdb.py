#!/usr/bin/env python3
import argparse
import os
import json
import logging
import typing as ty

from tqdm import tqdm

from retrieval import parse_antismash_json


# function to set up logger with stream handlers for stdout and to log file, you have to give path to log file
def setup_logger(log_file: str) -> None:
    """Set up logger."""
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

    # log to file
    file_handler = logging.FileHandler(log_file)
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    # log to stdout
    stream_handler = logging.StreamHandler()
    stream_handler.setLevel(logging.INFO)
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    return logger


def cli() -> argparse.Namespace:
    """Command line interface."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", type=str, required=True, help="path to dir with input jsons")
    parser.add_argument("--output", "-o", type=str, required=True, help="path to output dir")
    return parser.parse_args()


def parse_asdb_jsons(data: ty.Dict[str, ty.Any]) -> ty.Generator[ty.Any, None, None]:
    """Parse ASDB JSON files."""
    queries = parse_antismash_json(data, predict_specificities=False)
    for query_idx, query in enumerate(queries): 
        yield query_idx + 1, query


def main() -> None:
    "Main function."
    args = cli()
    logger = setup_logger("parse_asdb.log")
    
    # log input args
    logger.info(f"Input: {args.input}")

    # check if input and output are dirs
    if not os.path.isdir(args.input):
        print("Input is not a directory.")
        exit(1)

    if not os.path.isdir(args.output):
        print("Output is not a directory.")
        exit(1)

    # parse the input files, creating a geneator, and write the output
    file_paths = os.listdir(args.input)
    file_paths = sorted(file_paths)
    logger.info(f"Found {len(file_paths)} files.")

    logger.info("Parsing ASDB JSON files...")
    for file_path in tqdm(file_paths):
        try:
            if not file_path.endswith(".json"): continue
            with open(os.path.join(args.input, file_path), "r") as file: 
                data = json.load(file)
            for query_idx, query in parse_asdb_jsons(data):
                query["input_file"] = file_path
                query["protocluster"] = query_idx
                out_name = f"{query['record_id']}_protocluster{query_idx}.json"
                with open(os.path.join(args.output, out_name), "w") as out_file:
                    json.dump(query, out_file, indent=4)
        except Exception as e:
            logger.error(f"Error parsing {file_path}: {e}")
            continue
    
    logger.info("Done.")


if __name__ == "__main__":
    main()
