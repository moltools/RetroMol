#!/usr/bin/env python3
import argparse 


from .populate_database import connect_to_database


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, help="outdir deorphanization results")
    return parser.parse_args()


def main():
    args = cli()
    cur, conn = connect_to_database()

    # TODO
    # retrieve primary sequences eligible for deorphanization directly from database

    # TODO
    # for deorphanization:
    # filter on primary sequences related to protoclusters, not compounds
    # filter on organism (at least type and genus)
    # filter on length (somewhat, the iteratives are not interesting anyway -- prioritize longer and exact matches)
    # write out scored results with top 50 matches to file (include alignment)


if __name__ == "__main__":
    main()
