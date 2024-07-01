# -*- coding: utf-8 -*-

"""Command line interface for :mod:`retromol.retrosynthesis`."""

import argparse
import json
import logging
import multiprocessing as mp
import os
import typing as ty

from tqdm import tqdm

from retromol.retrosynthesis.chem import MolecularPattern, Molecule, ReactionRule
from retromol.retrosynthesis.helpers import timeout
from retromol.retrosynthesis.parsing import (
    Result,
    parse_mol,
    parse_molecular_patterns,
    parse_reaction_rules,
)

__all__ = ["main"]


def add_subparsers(parser: argparse._SubParsersAction) -> None:
    """Parse command line arguments.

    :param parser: The parser to add subparsers to.
    :type parser: argparse._SubParsersAction
    """
    subparser_single_mode = parser.add_parser("retrosynthesis_single")
    subparser_batch_mode = parser.add_parser("retrosynthesis_batch")

    repository_path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
    # fixtures_path = os.path.join(repository_path, "tests/fixtures")
    fixtures_path = fixtures_path = os.path.join("app", "src", "server", "data")

    subparser_single_mode.add_argument(
        "-r",
        "--reactions",
        type=str,
        required=False,
        default=os.path.join(fixtures_path, "reactions.json"),
        help="Path to JSON file containing reaction.",
    )
    subparser_single_mode.add_argument(
        "-m",
        "--monomers",
        type=str,
        required=False,
        default=os.path.join(fixtures_path, "monomers.json"),
        help="Path to JSON file containing monomer patterns.",
    )

    subparser_batch_mode.add_argument(
        "-r",
        "--reactions",
        type=str,
        required=False,
        default=os.path.join(fixtures_path, "reactions.json"),
        help="Path to JSON file containing reaction.",
    )
    subparser_batch_mode.add_argument(
        "-m",
        "--monomers",
        type=str,
        required=False,
        default=os.path.join(fixtures_path, "monomers.json"),
        help="Path to JSON file containing monomer patterns.",
    )

    subparser_single_mode.add_argument(
        "-i", "--input", type=str, required=True, help="Input SMILES (between quotes)."
    )
    subparser_single_mode.add_argument(
        "-o", "--output", type=str, required=True, help="Path to new JSON output file."
    )

    subparser_batch_mode.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="Path to input CSV or TSV file as 'name,SMILES' per line.",
    )
    subparser_batch_mode.add_argument(
        "-o",
        "--output",
        type=str,
        required=True,
        help="Path to new directory for JSON output files.",
    )
    subparser_batch_mode.add_argument(
        "-n",
        "--nproc",
        type=int,
        default=mp.cpu_count(),
        help="Number of processes to use.",
    )
    subparser_batch_mode.add_argument(
        "--header", action="store_true", help="Input file contains a header line."
    )


Record = ty.Tuple[Molecule, ty.List[ReactionRule], ty.List[MolecularPattern]]


@timeout(10)
def parse_mol_timed(record: Record) -> Result:
    """Parse a molecule with a timeout.

    :param record: The record containing the molecule, reaction rules, and
        molecular patterns.
    :type record: Record
    :return: The result of the parsing.
    :rtype: Result
    """
    logger = logging.getLogger(__name__)
    try:
        mol, reactions, monomers = record
        return parse_mol(mol, reactions, monomers)
    except Exception as e:
        logger.error(f"Failed to parse {mol.name} and raised {e.__class__.__name__}: {e}")
        return Result(mol.name, mol.compiled, False)


def main(args) -> None:
    """Driver function."""
    logger = logging.getLogger(__name__)

    # Parse reaction rules and molecular patterns.
    reactions_src = json.load(open(args.reactions, "r", encoding="utf-8"))
    reactions = parse_reaction_rules(reactions_src)

    monomers_src = json.load(open(args.monomers, "r", encoding="utf-8"))
    monomers = parse_molecular_patterns(monomers_src)

    logger.debug(f"Parsed {len(reactions)} reactions and {len(monomers)} monomers.")

    # Parse molecule in single mode.
    if args.mode == "retrosynthesis_single":
        mol = Molecule("input", args.input)
        record = (mol, reactions, monomers)
        result = parse_mol_timed(record)
        with open(args.output, "w", encoding="utf-8") as fo:
            fo.write(result.to_json())
        logger.debug("Processed single molecule.")
        exit(0)

    # Parse a batch of molecules.
    elif args.mode == "retrosynthesis_batch":

        # Parse input file.
        records = []
        with open(args.input, "r", encoding="utf-8") as fo:
            # Check file extension.
            sep = {"csv": ",", "tsv": "\t"}.get(args.input.split(".")[-1])
            if args.header:
                next(fo)
            for line in tqdm(fo, desc="Reading input file"):
                name, smiles = line.strip().split(sep)

                # Check if output dir already contains JSON file for entry.
                out_path = os.path.join(args.output, f"{name}.json")

                # Skip if JSON file already exists.
                if os.path.exists(out_path):
                    continue

                mol = Molecule(name, smiles)
                records.append((mol, reactions, monomers))

        logger.info(f"Added {len(records)} records to queue.")

        # Parse molecules in parallel.
        nproc = min(args.nproc, mp.cpu_count())
        logger.info(f"\nUsing {nproc} processes...")
        with mp.Pool(processes=nproc) as pool:
            for result in tqdm(pool.imap_unordered(parse_mol_timed, records)):
                out_path = os.path.join(args.output, f"{result.identifier}.json")
                with open(out_path, "w", encoding="utf-8") as fo:
                    fo.write(result.to_json())

        logger.info(f"Processed {len(records)} records.")
        exit(0)

    # Invalid mode.
    else:
        logger.error(f"Invalid mode: {args.mode}")
        exit(1)
