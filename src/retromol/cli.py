# -*- coding: utf-8 -*-

"""Command linwe interface for :mod:`retromol`."""

import argparse
import logging

import retromol.npkg.cli
import retromol.retrosynthesis.cli
from retromol.version import get_version

__all__ = ["main"]


def cli() -> argparse.Namespace:
    """Parse command line arguments.

    :return: The parsed command line arguments.
    :rtype: argparse.Namespace
    """
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {get_version()}",
        help="Show program's version number and exit.",
    )  # noqa: E501
    parser.add_argument(
        "-h",
        "--help",
        action="help",
        default=argparse.SUPPRESS,
        help="Show this help message and exit.",
    )  # noqa: E501
    parser.add_argument(
        "-l", "--logger-level", type=str, default="INFO", help="Logger verbosity level."
    )  # noqa: E501
    subparser = parser.add_subparsers(dest="mode", required=True)

    retromol.retrosynthesis.cli.add_subparsers(subparser)
    retromol.npkg.cli.add_subparsers(subparser)

    return parser.parse_args()


def main() -> None:
    """Driver function."""
    # Parse the command line arguments.
    args = cli()

    # Set up the logger.
    logging.basicConfig(
        level=args.logger_level, format="%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    # Dispatch to the appropriate mode.
    if args.mode in ["retrosynthesis_single", "retrosynthesis_batch"]:
        retromol.retrosynthesis.cli.main(args)

    elif args.mode in ["npkg_create", "npkg_purge", "npkg_add_protoclusters", "npkg_purge_protoclusters"]:
        retromol.npkg.cli.main(args)

    else:
        msg = f"Unknown mode: {args.mode}"
        logging.error(msg)
        raise ValueError(msg)

    exit(0)


if __name__ == "__main__":
    main()
