"""This module contains the command line interface for RetroMol."""

import argparse
import json
import logging
import os
import re
from collections import Counter
from contextlib import ExitStack
from datetime import datetime
from typing import Any

from tqdm import tqdm
from rdkit import RDLogger

from retromol.utils.logging import setup_logging, add_file_handler
from retromol.model.rules import RuleSet
from retromol.model.result import Result
from retromol.model.assembly_sampling import AssemblyConstraintError, AssemblyLimitError
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol_with_timeout
from retromol.io.streaming import run_retromol_stream, stream_sdf_records, stream_table_rows, stream_json_records
from retromol.io.json import dumps_compact, open_json_file
from retromol.visualization.reaction_graph import visualize_reaction_graph
from retromol.visualization.assembly_layout import MoleculeLayout


log = logging.getLogger(__name__)


RDLogger.DisableLog('rdApp.*')  # disable RDKit warnings


from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("retromol")
except PackageNotFoundError:
    # When running in a source checkout and haven’t installed it yet,
    # importlib.metadata might not find “retromol” in site‐packages.
    # Fallback to a hard-coded default or raise an error:
    __version__ = "0.0.0"


def add_rule_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--rxn-rules", type=str, default=None, help="path to reaction rules YAML (default: None, default ruleset)")
    parser.add_argument("--mxn-rules", type=str, default=None, help="path to matching rules YAML (default: None, default ruleset)")


def nonnegative_int(value: str) -> int:
    """Validate counts at argument parsing time, before running chemistry."""
    try:
        count = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a nonnegative integer") from None
    if count < 0:
        raise argparse.ArgumentTypeError("must be a nonnegative integer")
    return count


def coverage_fraction(value: str) -> float:
    """Validate an inclusive coverage threshold before running chemistry."""
    try:
        coverage = float(value)
    except ValueError:
        raise argparse.ArgumentTypeError("must be a finite fraction between 0 and 1") from None
    if not 0 <= coverage <= 1:
        raise argparse.ArgumentTypeError("must be a finite fraction between 0 and 1")
    return coverage


def cli() -> argparse.Namespace:
    """
    Parse command line arguments.

    :return: Parsed command line arguments.
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
    single_parser.add_argument(
        "--max-assembly-graphs", type=nonnegative_int, default=1, metavar="N",
        help="maximum assembly drawings: 0 skips, 1 prefers the default assembly if eligible, "
             "larger values sample distinct eligible assemblies (default: 1)",
    )
    single_parser.add_argument(
        "--min-coverage", type=coverage_fraction, default=None, metavar="FRACTION",
        help="minimum identified heavy-atom coverage from 0 to 1; "
             "default: only assemblies with the highest attainable coverage; 0 allows all views",
    )
    single_parser.add_argument(
        "--seed", type=int, default=42,
        help="assembly sampling seed (default: 42)",
    )
    single_parser.add_argument(
        "--compression", choices=["none", "gzip"], default="none",
        help="result JSON compression (default: none; gzip writes result.json.gz)",
    )

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
    batch_parser.add_argument("--jsonl-path", type=str, default=None, help="path to results JSONL (.gz appended for gzip; default: <outdir>/results.jsonl.gz)")
    batch_parser.add_argument(
        "--compression", choices=["none", "gzip"], default="gzip",
        help="result compression (default: gzip; use none for plain JSONL/JSON)",
    )
    batch_parser.add_argument("--no-tqdm", action="store_true", help="disable progress bars for lowest overhead")
    batch_parser.add_argument("--rdkit-fast", action="store_true", help="use fast SDF parse (sanitize=False, removeHs=True); we'll sanitize only when needed")

    # Only read when input type is table
    batch_parser.add_argument("--separator", type=str, choices=["comma", "tab"], default="comma", help="separator for table file (default: ',')")
    batch_parser.add_argument("--id-col", type=str, default="inchikey", help="name of the column containing InChIKeys (default: 'inchikey')")
    batch_parser.add_argument("--smiles-col", type=str, default="smiles", help="name of the column containing SMILES strings (default: 'smiles')")

    # Batch mode also allows for parallel processing
    batch_parser.add_argument("-w", "--workers", type=int, default=1, help="number of worker processes to use (default: 1)")

    # Allow custom rule sets for both single and batch parsers
    add_rule_args(single_parser)
    add_rule_args(batch_parser)

    args = parser.parse_args()
    if (args.mode == "batch" and args.jsonl_path and args.jsonl_path.lower().endswith(".gz")
            and args.compression == "none"):
        parser.error("--jsonl-path ending in .gz requires --compression gzip")
    return args


def _open_jsonl(outdir: str, jsonl_path: str | None, compression: str = "gzip") -> tuple[Any, str]:
    """
    Open a plain or gzip JSONL stream for appending results.
    
    :param outdir: str: output directory 
    :param jsonl_path: str | None: path to JSONL file, or None to use default
    :param compression: gzip (default) or none; append .gz when needed
    :return: tuple[file handle, path]: opened file handle and the path used
    """
    path = jsonl_path or os.path.join(outdir, "results.jsonl")
    if compression == "gzip" and not path.lower().endswith(".gz"):
        path += ".gz"
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    # One continuous stream per run, not one gzip member per compound. This
    # lets compression reuse repeated structures/metadata across records.
    return open_json_file(path, "at"), path


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
    ruleset = RuleSet.load_from_files(
        reaction_rules_path=args.rxn_rules,
        matching_rules_path=args.mxn_rules,
        match_stereochemistry=args.c
    )
    log.info(f"loaded default ruleset: {ruleset}")

    result_counts: Counter[str] = Counter()

    # Single mode
    if args.mode == "single":
        submission = Submission(args.smiles, props={})
        result: Result = run_retromol_with_timeout(submission, ruleset)
        log.info(f"result: {result}")

        # Write out result to file and then read back in again for visualization (test I/O)
        result_dict = result.to_dict()
        result_path = os.path.join(args.outdir, "result.json" + (".gz" if args.compression == "gzip" else ""))
        with open_json_file(result_path, "wt") as f:
            f.write(dumps_compact(result_dict))

        with open_json_file(result_path, "rt") as f:
            result_data = json.load(f)
        result2 = Result.from_dict(result_data)

        if args.max_assembly_graphs == 0:
            assemblies = []
            log.info("Skipping assembly drawings (--max-assembly-graphs 0)")
        else:
            space = result2.assembly_space()
            try:
                threshold = args.min_coverage
                default_coverage = result2.calculate_coverage() if args.max_assembly_graphs == 1 else None
                if threshold is None:
                    # A fully covered default already proves the maximum and
                    # avoids enumeration for the usual single-drawing case.
                    threshold = 1.0 if default_coverage == 1.0 else space.max_coverage
                if default_coverage is not None and default_coverage >= threshold:
                    assemblies = [result2.assembly_graph]
                else:
                    assemblies = space.sample(args.max_assembly_graphs, seed=args.seed, min_coverage=threshold)
                    log.info(
                        "Drawing %d of %d eligible assemblies (%d total; minimum coverage %.2f%%; seed=%d)",
                        len(assemblies), space.count(min_coverage=threshold), space.count(), 100 * threshold, args.seed,
                    )
                    if not assemblies:
                        if not space.frontiers:
                            log.warning("No explored assemblies preserve the identified terminal units; skipping drawings.")
                        else:
                            log.warning("No assemblies meet minimum coverage %.2f%%; highest attainable coverage is %.2f%%.",
                                        100 * threshold, 100 * space.max_coverage)
            except (AssemblyLimitError, AssemblyConstraintError) as exc:
                log.error("Cannot sample assembly drawings: %s. The reaction graph is saved in %s.", exc, result_path)
                raise SystemExit(1) from None

        # Rerunning into the same directory must not leave earlier, now-invalid
        # projections visible beside the new drawings. Only our generated names
        # are replaced; user-named images in the directory are left alone.
        for entry in os.scandir(args.outdir):
            if entry.is_file() and re.fullmatch(r"assembly_graph(?:_[0-9]{3,})?\.png", entry.name):
                os.remove(entry.path)

        drawing_layout = MoleculeLayout.from_mol(result2.root.mol) if assemblies else None
        for index, assembly in enumerate(assemblies, start=1):
            filename = "assembly_graph.png" if args.max_assembly_graphs == 1 else f"assembly_graph_{index:03d}.png"
            assembly.draw(show_unassigned=True, savepath=os.path.join(args.outdir, filename), layout=drawing_layout)
            log.info("%s coverage: %.2f%%", filename, 100 * result2.calculate_coverage(assembly))

        # Visualize reaction graph
        visualize_reaction_graph(
            result2.reaction_graph,
            html_path=os.path.join(args.outdir, "reaction_graph.html"),
            root_enc=result2.root_enc
        )

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

        result_counts = Counter()
        processed_in_current_batch = 0
        # Close streams even if parsing/writing is interrupted, so completed
        # records and the gzip footer are flushed to disk.
        with ExitStack() as resources:
            pbar_outer = resources.enter_context(tqdm(desc="Batches", unit="batch", disable=args.no_tqdm))
            pbar_inner = resources.enter_context(tqdm(desc="Processed", unit="mol", disable=args.no_tqdm))
            jsonl_fh = None
            if args.results == "jsonl":
                jsonl_fh, jsonl_path = _open_jsonl(args.outdir, args.jsonl_path, args.compression)
                resources.enter_context(jsonl_fh)
                log.info(f"Appending results to JSONL file at: {jsonl_path}")

            for record_index, evt in enumerate(run_retromol_stream(
                ruleset=ruleset,
                row_iter=source_iter,
                smiles_col=smiles_col,
                workers=args.workers,
                batch_size=args.batch_size,
                pool_chunksize=args.pool_chunksize,
                maxtasksperchild=args.maxtasksperchild,
            ), start=1):
                if evt.error is not None:
                    log.error(evt.error)
                    result_counts["errors"] += 1
                elif evt.result is not None:
                    serialized = dumps_compact(evt.result)
                    if jsonl_fh is not None:
                        jsonl_fh.write(serialized + "\n")
                    else:
                        # Numbered files preserve repeated input IDs. Exclusive
                        # creation avoids silently replacing a previous run.
                        suffix = ".json.gz" if args.compression == "gzip" else ".json"
                        path = os.path.join(args.outdir, f"result_{record_index:09d}{suffix}")
                        with open_json_file(path, "xt") as f:
                            f.write(serialized)
                    result_counts["successes"] += 1
                else:
                    log.error("received empty result with no error message")
                    result_counts["failures"] += 1

                pbar_inner.update(1)
                processed_in_current_batch += 1
                if processed_in_current_batch >= args.batch_size:
                    pbar_outer.update(1)
                    processed_in_current_batch = 0

            if processed_in_current_batch > 0:
                pbar_outer.update(1)

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
