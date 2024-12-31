#!/usr/bin/env python
import os
import psycopg2
import argparse 
from modules.step1_add_npatlas import add_npatlas
from modules.step2_add_donphan import add_donphan
from modules.step3_add_retromol_results import add_retromol_results


def cli():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--npatlas", required=False, help="path to the NPAtlas database json file")
    parser.add_argument("--donphan", required=False, help="path to the DONPHAN database csv file")
    parser.add_argument("--retromol", required=False, help="path to RetroMol results folder for NPAtlas compounds")
    return parser.parse_args()


def connect_to_database():
    """Connect to the database and return the connection and cursor."""
    conn = psycopg2.connect(
        dbname=os.getenv("DB_NAME", "retromol"),
        user=os.getenv("DB_USER", "davidmeijer"),
        password=os.getenv("DB_PASSWORD", "postgres"),
        host=os.getenv("DB_HOST", "localhost"),
        port=os.getenv("DB_PORT", "5432")
    )
    return conn.cursor(), conn


def main():
    args = cli()
    path_npatlas = args.npatlas
    path_donphan = args.donphan
    path_retromol = args.retromol
    path_mibig_database = args.mibig_database
    path_mibig_jsons = args.mibig_jsons

    cur, conn = connect_to_database()

    # add the NPAtlas database to the database
    if path_npatlas:
        add_npatlas(cur, path_npatlas)
        conn.commit()
    
    # add the DONPHAN database to the database
    if path_donphan:
        add_donphan(cur, path_donphan)
        conn.commit()

    # add the RetroMol results to the database
    if path_retromol:
        # NOTE: running this module twice, will add the same primary sequences twice
        add_retromol_results(cur, path_retromol)
        conn.commit()

    # count the number of compound records in the database
    cur.execute("SELECT COUNT(*) FROM compounds")
    count = cur.fetchone()[0]
    print(f"number of compound records in the database: {count}")

    # count the number of organism records in the database
    cur.execute("SELECT COUNT(*) FROM organisms")
    count = cur.fetchone()[0]
    print(f"number of organism records in the database: {count}")

    # count the number of biosynthetic class records in the database
    cur.execute("SELECT COUNT(*) FROM biosynthetic_classes")
    count = cur.fetchone()[0]
    print(f"number of biosynthetic_class records in the database: {count}")

    # print top 3 most associated bioactivities
    cur.execute("SELECT cb.bioactivity_type, COUNT(DISTINCT cb.compound_id) AS compound_count FROM compounds_bioactivities AS cb GROUP BY cb.bioactivity_type ORDER BY compound_count DESC LIMIT 3;")
    for bioactivity, count in cur.fetchall():
        print(f"bioactivity: {bioactivity}, count: {count}")

    # count the number of primary_sequence records in the database
    cur.execute("SELECT COUNT(*) FROM primary_sequences")
    count = cur.fetchone()[0]
    print(f"number of primary_sequence records in the database: {count}")

    cur.close()
    conn.close()


if __name__ == "__main__":
    main()
