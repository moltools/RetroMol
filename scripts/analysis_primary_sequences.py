#!/usr/bin/env python3
import argparse 


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", required=True, type=str, help="path to retromol results folder")
    parser.add_argument("-a", required=True, type=str, help="path to original retromol input file with annotations")
    parser.add_argument("-o", required=True, type=str, help="path to output folder")
    return parser.parse_args()


def main() -> None:
    args = cli()


if __name__ == "__main__":
    main()
