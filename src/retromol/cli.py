"""This module contains the command line interface for RetroMol."""

import argparse
import json
import logging
import os
import time
from collections import Counter
from datetime import datetime
from typing import Any

from tqdm import tqdm
from rdkit import RDLogger

from retromol.version import __version__
from retromol.utils.logging import setup_logging, add_file_handler
from retromol.model.rules import RuleSet
from retromol.model.result import Result
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol_with_timeout
from retromol.io.streaming import run_retromol_stream, stream_sdf_records, stream_table_rows, stream_json_records
from retromol.chem.mol import encode_mol
from retromol.visualization.reaction_graph import visualize_reaction_graph
from retromol.fingerprint.fingerprint import FingerprintGenerator


log = logging.getLogger(__name__)


RDLogger.DisableLog('rdApp.*')  # disable RDKit warnings


def cli() -> argparse.Namespace:
    """
    Parse command line arguments.

    :return: parsed command line arguments
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("-o", "--outdir", type=str, required=True, help="output directory for results")
    parser.add_argument("-l", "--log-level", type=str, choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="INFO", help="logging level (default: INFO)")

    parser.add_argument("-h", "--help", action="help", help="show cli options")
    parser.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    parser.add_argument("-c", action="store_true", help="match stereochemistry in the input SMILES (default: False)")

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
    batch_parser.add_argument("--batch-size", type=int, default=2000, help="max tasks buffered before dispatch (default: 2000)")
    batch_parser.add_argument("--chunksize", type=int, default=20000, help="rows per CSV/TSV chunk (default: 20000)")
    batch_parser.add_argument("--pool-chunksize", type=int, default=50, help="chunksize hint for imap_unordered (default: 50)")
    batch_parser.add_argument("--maxtasksperchild", type=int, default=2000, help="recycle worker after N tasks (default: 2000)")
    batch_parser.add_argument("--results", choices=["files", "jsonl"], default="jsonl", help="write each result to a file or append to JSONL (default: jsonl)")
    batch_parser.add_argument("--jsonl-path", type=str, default=None, help="path to results jsonl (default: <outdir>/results.jsonl)")
    batch_parser.add_argument("--no-tqdm", action="store_true", help="disable progress bars for lowest overhead")
    batch_parser.add_argument("--rdkit-fast", action="store_true", help="use fast SDF parse (sanitize=False, removeHs=True); we'll sanitize only when needed")

    # Only read when input type is table
    batch_parser.add_argument("--separator", type=str, choices=["comma", "tab"], default="comma", help="separator for table file (default: ',')")
    batch_parser.add_argument("--id-col", type=str, default="inchikey", help="name of the column containing InChIKeys (default: 'inchikey')")
    batch_parser.add_argument("--smiles-col", type=str, default="smiles", help="name of the column containing SMILES strings (default: 'smiles')")

    # Batch mode also allows for parallel processing
    batch_parser.add_argument("-w", "--workers", type=int, default=1, help="number of worker processes to use (default: 1)")

    return parser.parse_args()


def _open_jsonl(outdir: str, jsonl_path: str | None) -> tuple[Any, str]:
    """
    Open a JSONL file for appending results.
    
    :param outdir: str: output directory 
    :param jsonl_path: str | None: path to JSONL file, or None to use default
    :return: tuple[file handle, path]: opened file handle and the path used
    """
    path = jsonl_path or os.path.join(outdir, "results.jsonl")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    return open(path, "a", buffering=1), path  # line-buffered


def main() -> None:
    """
    Main entry point for the CLI.
    """
    start_time = datetime.now()

    # Parse command line arguments and set up logging
    args = cli()

    # Create output directory if it doesn't exist
    os.makedirs(args.outdir, exist_ok=True)

    # Setup logging
    setup_logging(level=args.log_level)

    # If log file exists, remove it
    log_fp = os.path.join(args.outdir, "retromol.log")
    if os.path.exists(log_fp):
        os.remove(log_fp)
    
    # Add file handler to log to file
    add_file_handler(log_fp, level=args.log_level)
    
    # Log command line arguments
    log.info("command line arguments:")
    for arg, val in vars(args).items():
        log.info(f"\t{arg}: {val}")

    # Load default ruleset
    ruleset = RuleSet.load_default(match_stereochemistry=args.c)
    log.info(f"loaded default ruleset: {ruleset}")

    result_counts: Counter[str] = Counter()

    # Single mode
    if args.mode == "single":
        submission = Submission(args.smiles, props={})
        result: Result = run_retromol_with_timeout(submission, ruleset)
        log.info(f"result: {result}")

        # Write out result to file and then read back in again for visualization (test I/O)
        result_dict = result.to_dict()
        with open(os.path.join(args.outdir, "result.json"), "w") as f:
            json.dump(result_dict, f, indent=4)

        with open(os.path.join(args.outdir, "result.json"), "r") as f:
            result_data = json.load(f)
        result2 = Result.from_dict(result_data)

        # Report on coverage as percentage of tags identified
        coverage = result2.calculate_coverage()
        log.info(f"coverage: {coverage:.2%}")

        # Get linear readout; draw assembly graph
        linear_readout = result2.linear_readout
        out_assembly_graph_fig = os.path.join(args.outdir, "assembly_graph.png")
        linear_readout.assembly_graph.draw(show_unassigned=True, savepath=out_assembly_graph_fig)
        log.info(f"linear readout: {linear_readout}")

        # Visualize reaction graph
        root = encode_mol(result2.submission.mol)
        visualize_reaction_graph(
            result2.reaction_graph,
            html_path=os.path.join(args.outdir, "reaction_graph.html"),
            root_enc=root
        )

        # Calculate a fingerprint for the molecule based on matched rules
        t0 = time.time()
        ruleset = RuleSet.load_default()
        matching_rules = ruleset.matching_rules
        generator = FingerprintGenerator(matching_rules)
        fp = generator.fingerprint_from_result(result, num_bits=512, counted=True)
        t1 = time.time()
        log.info(f"generated fingerprint with {len(fp)} bits set in {t1 - t0:.2f} seconds")
        log.info(f"fingerprint contains {sum(1 for v in fp if v != 0)} non-zero items")

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

        # Results saved into JSONL format to limit file operations
        jsonl_fh = None
        jsonl_path = None
        if args.results == "jsonl":
            jsonl_fh, jsonl_path = _open_jsonl(args.outdir, args.jsonl_path)
            log.info(f"Appending results to JSONL file at: {jsonl_path}")

        result_counts = Counter()

        processed_in_current_batch = 0

        for evt in run_retromol_stream(
            ruleset=ruleset,
            row_iter=source_iter,
            smiles_col=smiles_col,
            workers=args.workers,
            batch_size=args.batch_size,
            pool_chunksize=args.pool_chunksize,
            maxtasksperchild=args.maxtasksperchild,
        ):
            # evt has: result (dict or None) and error (str or None)
            if evt.error is not None:
                log.error(evt.error)
                result_counts["errors"] += 1
            elif evt.result is not None:
                # Result is already serialized as dict
                jsonl_fh.write(json.dumps(evt.result) + "\n")
                result_counts["successes"] += 1
            else:
                log.error("received empty result with no error message")
                result_counts["failures"] += 1

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

        log.info(f"Streaming complete. Summary: {dict(result_counts)}")

    else:
        log.error("either --smiles or --database must be provided")

    log.info("processing complete")
    log.info(f"summary of results: {dict(result_counts)}")

    # Wrap up
    end_time = datetime.now()
    run_time = end_time - start_time
    log.info(f"start time: {start_time}, end time: {end_time}, run time: {run_time}")
    log.info("goodbye")


if __name__ == "__main__":
    main()
