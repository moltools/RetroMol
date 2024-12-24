#!/usr/bin/env python
import os
import psycopg2
import argparse 
from modules.step1_add_npatlas import add_npatlas
from modules.step2_add_donphan import add_donphan
from modules.step3_add_retromol_results import add_retromol_results
from database.modules.step4_add_mibig_depr import add_mibig


def cli():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--npatlas", required=False, help="path to the NPAtlas database json file")
    parser.add_argument("--donphan", required=False, help="path to the DONPHAN database csv file")
    parser.add_argument("--retromol", required=False, help="path to RetroMol results folder for NPAtlas compounds")
    parser.add_argument("--mibig-database", required=False, help="path to the MIBiG database json files folder (download from MIBiG website)")
    parser.add_argument("--mibig-jsons", required=False, help="path to the MIBiG jsons folder (download with 'get_mibig_jsons.py' script)")
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

    # add the MIBiG database to the database
    if path_mibig_database and path_mibig_jsons:
        add_mibig(cur, path_mibig_database, path_mibig_jsons)
        conn.commit()
    elif path_mibig_database or path_mibig_jsons:
        raise ValueError("both mibig database and mibig jsons must be provided")
    
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

    # count the number of protocluster records in the database
    cur.execute("SELECT COUNT(*) FROM protoclusters")
    count = cur.fetchone()[0]
    print(f"number of protocluster records in the database: {count}")

    # count the number of compound records linked to primary_sequences
    cur.execute("SELECT COUNT(DISTINCT primary_sequence_id) AS total_linked_sequences FROM compounds_primary_sequences;")
    count = cur.fetchone()[0]
    print(f"number of primary_sequences linked to compounds: {count}")

    # count the number of protocluster records linked to primary_sequences
    cur.execute("SELECT COUNT(DISTINCT primary_sequence_id) AS total_linked_sequences FROM protoclusters_primary_sequences;")
    count = cur.fetchone()[0]
    print(f"number of primary_sequences linked to protoclusters: {count}")

    cur.close()
    conn.close()


if __name__ == "__main__":
    main()
