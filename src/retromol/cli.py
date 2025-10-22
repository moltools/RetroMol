"""This module contains the command line interface for RetroMol."""

import argparse
import json
import logging
import os
import os.path as osp
from collections import Counter
from datetime import datetime
from typing import Any

import yaml
from tqdm import tqdm

from retromol.api import run_retromol_with_timeout
from retromol.config import LOGGER_LEVEL, LOGGER_NAME
from retromol.drawing import draw_result
from retromol.io import Input as RetroMolInput
from retromol.io import Result
from retromol.readout import linear_readout_with_timeout
from retromol.rules import (
    get_path_default_matching_rules,
    get_path_default_reaction_rules,
    get_path_default_wave_config,
    load_rules_from_files,
)
from retromol.streaming import run_retromol_stream, stream_json_records, stream_sdf_records, stream_table_rows
from retromol.version import __version__

logger = logging.getLogger(LOGGER_NAME)
logger.setLevel(LOGGER_LEVEL)


def cli() -> argparse.Namespace:
    """
    Parse command line arguments.

    :return: parsed command line arguments
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-o", "--outdir", type=str, required=True, help="output directory for results")

    parser.add_argument("-h", "--help", action="help", help="show cli options")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")

    parser.add_argument(
        "-rx",
        "--reaction-rules",
        type=str,
        required=False,
        default=get_path_default_reaction_rules(),
        help="path to reaction rules yaml",
    )
    parser.add_argument(
        "-rm",
        "--matching-rules",
        type=str,
        required=False,
        default=get_path_default_matching_rules(),
        help="path to matching rules yaml",
    )
    parser.add_argument(
        "-wc",
        "--wave-config",
        type=str,
        required=False,
        default=get_path_default_wave_config(),
        help="path to wave configuration yaml",
    )

    # Flags
    parser.add_argument(
        "-C",
        "--matchstereochem",
        action="store_true",
        help="match stereochemistry in the input SMILES (default: False)",
    )
    parser.add_argument("-V", "--verbose", action="store_true", help="enable verbose output")
    parser.add_argument(
        "-D",
        "--check-duplicates",
        action="store_true",
        help="check for duplicate items in matching rules",
    )

    # Create two subparsers 'single' and 'batch'
    subparsers = parser.add_subparsers(dest="mode", required=True)
    single_parser = subparsers.add_parser("single", help="process a single compound")
    batch_parser = subparsers.add_parser("batch", help="process a batch of compounds")

    # For 'single' mode user should just give a SMILES as input
    single_parser.add_argument("-s", "--smiles", type=str, help="SMILES string of the compound to process")

    # For 'batch' mode user should provide a path to an SDF file
    input_group = batch_parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-s", "--sdf", type=str, help="path to an SDF file containing compounds to process")
    input_group.add_argument("-t", "--table", type=str, help="path to a CSV/TSV file containing compounds to process")
    input_group.add_argument("-j", "--json", type=str, help="path to a JSONL file containing compounds to process")
    batch_parser.add_argument(
        "--batch-size",
        type=int,
        default=2000,
        help="max tasks buffered before dispatch (default: 2000)",
    )
    batch_parser.add_argument("--chunksize", type=int, default=20000, help="rows per CSV/TSV chunk (default: 20000)")
    batch_parser.add_argument(
        "--pool-chunksize",
        type=int,
        default=50,
        help="chunksize hint for imap_unordered (default: 50)",
    )
    batch_parser.add_argument(
        "--maxtasksperchild",
        type=int,
        default=2000,
        help="recycle worker after N tasks (default: 2000)",
    )
    batch_parser.add_argument(
        "--results",
        choices=["files", "jsonl"],
        default="jsonl",
        help="write each result to a file or append to JSONL (default: jsonl)",
    )
    batch_parser.add_argument(
        "--jsonl-path",
        type=str,
        default=None,
        help="path to results jsonl (default: <outdir>/results.jsonl)",
    )
    batch_parser.add_argument("--no-tqdm", action="store_true", help="disable progress bars for lowest overhead")
    batch_parser.add_argument(
        "--rdkit-fast",
        action="store_true",
        help="use fast SDF parse (sanitize=False, removeHs=True); we’ll sanitize only when needed",
    )

    # Only read when input type is table
    batch_parser.add_argument(
        "--separator",
        type=str,
        choices=["comma", "tab"],
        default="comma",
        help="separator for table file (default: ',')",
    )
    batch_parser.add_argument(
        "--id-col",
        type=str,
        default="inchikey",
        help="name of the column containing InChIKeys (default: 'inchikey')",
    )
    batch_parser.add_argument(
        "--smiles-col",
        type=str,
        default="smiles",
        help="name of the column containing SMILES strings (default: 'smiles')",
    )

    # Batch mode also allows for parallel processing
    batch_parser.add_argument(
        "-w",
        "--workers",
        type=int,
        default=1,
        help="number of worker processes to use (default: 1)",
    )

    return parser.parse_args()


def setup_logger(log_file_path: str, verbose: bool) -> logging.Logger:
    """
    Sets up a logger that ONLY uses the handlers you attach,
    and won't bubble up to the root logger.

    :param log_file_path: path to the log file
    :param verbose: if True, also log to stdout
    :return: configured logger instance
    """
    # Remove old log file
    if os.path.exists(log_file_path):
        os.remove(log_file_path)

    # Remove any handlers that were already attached
    if logger.hasHandlers():
        logger.handlers.clear()

    # Common formatter
    fmt = logging.Formatter("[%(asctime)s][%(levelname)s]: %(message)s")

    # File handler
    fh = logging.FileHandler(log_file_path)
    fh.setLevel(LOGGER_LEVEL)
    fh.setFormatter(fmt)
    logger.addHandler(fh)

    # Only add stream handler if verbose is True
    if verbose:
        sh = logging.StreamHandler()
        sh.setLevel(LOGGER_LEVEL)
        sh.setFormatter(fmt)
        logger.addHandler(sh)

    return logger


def _open_jsonl(outdir: str, jsonl_path: str | None) -> tuple[Any, str]:
    path = jsonl_path or os.path.join(outdir, "results.jsonl")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return open(path, "a", buffering=1), path  # line-buffered


def _write_result_file(outdir: str, inchikey: str, payload: dict[str, Any] | None) -> None:
    with open(os.path.join(outdir, f"result_{inchikey}.json"), "w") as f:
        json.dump(payload, f, indent=0)  # indent=0 faster than 4


def main() -> None:
    """
    Main entry point for the CLI.
    """
    # Parse command line arguments and set up logging
    start_time = datetime.now()
    args = cli()
    os.makedirs(args.outdir, exist_ok=True)
    log_file_path = osp.join(args.outdir, "_retromol.log")  # add underscore to make log file appear at top of folder
    logger = setup_logger(log_file_path, args.verbose)
    logger.debug(f"command line arguments: {args}")

    # Load rules from files
    path_reaction_rules = args.reaction_rules
    path_matching_rules = args.matching_rules
    rule_set = load_rules_from_files(path_reaction_rules, path_matching_rules)
    logger.info(f"Loaded rule set: {rule_set}")

    # Check for duplicates if flag is set
    if args.check_duplicates:
        rule_set.check_for_duplicates()
        logger.info("Checked for duplicates in the rule set. Please remove duplicates for better performance.")

    # Load wave configuration
    with open(args.wave_config) as f:
        wave_configs = yaml.safe_load(f)
    logger.info(f"Loaded {len(wave_configs)} wave configuration(s) from {args.wave_config}")
    logger.debug(f"Loaded wave configuration: {wave_configs}")

    result_counts: Counter[str] = Counter()

    # Single mode
    if args.mode == "single":
        mol = RetroMolInput("target", args.smiles, props={})
        result: Result = run_retromol_with_timeout(mol, rule_set, wave_configs, args.matchstereochem)
        logger.info(f"Result: {result}")
        serialized_result = result.serialize()
        with open(osp.join(args.outdir, "result.json"), "w") as f:
            json.dump(serialized_result, f, indent=4)

        summary = result.summarize_by_depth()
        with open(osp.join(args.outdir, "summary.json"), "w") as f:
            json.dump(summary, f, indent=4)

        best_total_cov = result.best_total_coverage(round_to=2)
        logger.info(f"Best total coverage (rounded to 2 decimals): {best_total_cov}")

        draw_result(result, args.outdir, background_color="#fffaf6")

        logger.info("Writing linear readouts to log...")
        linear_readout = linear_readout_with_timeout(result)
        for level in linear_readout["levels"]:
            for path in level["strict_paths"]:
                ordered_monomers = path["ordered_monomers"]
                names = [m["identity"] for m in ordered_monomers]
                logger.info(f"{level['dfs_index']} - {level['depth']} - {names}")

        result_counts["successes"] += 1

    # Batch mode
    elif args.mode == "batch":
        id_col = args.id_col
        smiles_col = args.smiles_col
        separator = "," if args.separator == "comma" else "\t"

        # Choose source iterator (streamed, chunked)
        if args.sdf:
            source_iter = stream_sdf_records(args.sdf, fast=args.rdkit_fast)
        elif args.table:
            source_iter = stream_table_rows(args.table, sep=separator, chunksize=args.chunksize)
        else:
            source_iter = stream_json_records(args.json)

        # Progress bars: outer ~batches, inner = molecules processed
        pbar_outer = tqdm(desc="Batches", unit="batch", disable=args.no_tqdm)
        pbar_inner = tqdm(desc="Processed", unit="mol", disable=args.no_tqdm)

        # Result sink (same behavior as before)
        jsonl_fh = None
        jsonl_path = None
        if args.results == "jsonl":
            jsonl_fh, jsonl_path = _open_jsonl(args.outdir, args.jsonl_path)
            logger.info(f"Appending results to JSONL file at: {jsonl_path}")

        result_counts = Counter()

        processed_in_current_batch = 0

        for evt in run_retromol_stream(
            # Config
            rule_set=rule_set,  # already loaded above
            wave_configs=wave_configs,  # already loaded above
            match_stereo=args.matchstereochem,
            # Data & schema
            row_iter=source_iter,
            id_col=id_col,
            smiles_col=smiles_col,
            # Concurrency knobs
            workers=args.workers,
            batch_size=args.batch_size,
            pool_chunksize=args.pool_chunksize,
            maxtasksperchild=args.maxtasksperchild,
        ):
            # evt has: inchikey, result (dict or None), error (str or None)
            input_id = evt.inchikey
            if evt.error is not None:
                logger.error(f"Error {input_id}: {evt.error}")
                result_counts["errors"] += 1
            else:
                if args.results == "files":
                    _write_result_file(args.outdir, input_id, evt.result)
                else:
                    jsonl_fh.write(json.dumps({"inchikey": input_id, "result": evt.result}) + "\n")
                result_counts["successes"] += 1

            # Progress
            pbar_inner.update(1)
            processed_in_current_batch += 1
            if processed_in_current_batch >= args.batch_size:
                pbar_outer.update(1)
                processed_in_current_batch = 0

        # If there was a final partial batch, tick the outer bar once more
        if processed_in_current_batch > 0:
            pbar_outer.update(1)

        pbar_inner.close()
        pbar_outer.close()

        if jsonl_fh:
            jsonl_fh.close()

        logger.info(f"Streaming complete. Summary: {dict(result_counts)}")

    else:
        logger.error("Either --smiles or --database must be provided.")

    logger.info(f"Processing complete. Summary of results: {dict(result_counts)}")

    # Wrap up
    end_time = datetime.now()
    run_time = end_time - start_time
    logger.info(f"start time: {start_time}, end time: {end_time}, run time: {run_time}")
    logger.info("Goodbye.")


if __name__ == "__main__":
    main()
