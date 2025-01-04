#!/usr/bin/env python3
import argparse 
import os
import csv
import io
import json

from tqdm import tqdm

from populate_database import connect_to_database


# - MIBiG BGC y/n
# - best matching score to MIBiG BGC
# - MIBiG compound name(s) of strong match(es)


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", "-i", type=str, required=True, help="path to dir with deorphanization results")
    parser.add_argument("--mibig", "-m", type=str, required=True, help="path to MIBiG database")
    parser.add_argument("--output", "-o", type=str, required=True, help="path to output file")
    return parser.parse_args()


def get_compound_name(cur, compound_id):
    query = f"""
    SELECT compounds.name FROM compounds WHERE compounds.compound_id = '{compound_id}';
    """
    cur.execute(query)
    return cur.fetchone()[0]


def make_npatlas_url(compound_id):
    return f"https://www.npatlas.org/explore/compounds/{compound_id}"


def parse_mibig(path):
    bgc_to_compounds = {}
    file_paths = [os.path.join(path, f) for f in os.listdir(path) if f.endswith(".json")]
    for file_path in tqdm(file_paths):
        compound_names = []
        with open(file_path, "r") as f:
            data = json.load(f)
        accession = data["accession"]
        for compound in data["compounds"]:
            compound_names.append(compound["name"])
        bgc_to_compounds[accession] = "|".join(compound_names)
    return bgc_to_compounds


def main():
    args = cli()
    file_paths = [os.path.join(args.input, f) for f in os.listdir(args.input) if f.endswith(".csv")]
    print(f"Number of files: {len(file_paths)}")

    cur, conn = connect_to_database()

    out_file = open(args.output, "w")
    out_file.write("compound_id,compound_url,compound_name,highest_ranking_protocluster_score,length_compound_motif_code,length_protocluster_motif_code,highest_ranking_bgc_rank,highest_ranking_bgc,highest_ranking_bgc_score,highest_ranking_bgc_compound_names\n")

    high_scorers = 0
    records = []

    bgc_to_compounds = parse_mibig(args.mibig)
    
    for file_path in tqdm(file_paths):
        highest_ranking_score = 0.0
        length_compound_motif_code = 0
        length_protocluster_motif_code = 0

        highest_ranking_bgc_rank = ""
        highest_ranking_bgc = ""
        highest_ranking_bgc_compound_name = ""
        highest_ranking_bgc_score = 0.0

        with open(file_path, "r") as f:
            f.readline() # skip header
            for line_idx, line in enumerate(f):
                csv_file = io.StringIO(line)
                reader = csv.reader(csv_file)
                parsed_fields = next(reader)
                (
                    compound_id, 
                    protocluster_id, 
                    protocluster_input_file, 
                    score, 
                    _, 
                    _, 
                    _, 
                    _,
                    compound_motif_code,
                    protocluster_motif_code
                # ) = line.strip().split(",") # only need the first line
                ) = parsed_fields

                if line_idx == 0:
                    highest_ranking_score = float(score)
                    length_compound_motif_code = len(compound_motif_code.split("|"))
                    length_protocluster_motif_code = len(protocluster_motif_code.split("|"))

                if "BGC" in protocluster_input_file and len(highest_ranking_bgc) == 0:
                    highest_ranking_bgc_rank = str(line_idx + 1)
                    highest_ranking_bgc = protocluster_input_file.split(".", 1)[0]
                    highest_ranking_bgc_score = float(score)
                    highest_ranking_bgc_compound_name = bgc_to_compounds.get(highest_ranking_bgc, "Unknown")
            

        # only print 
        if highest_ranking_score >= 0.0:
            high_scorers += 1
            compound_name = get_compound_name(cur, compound_id)
            records.append((
                compound_id, 
                compound_name, 
                highest_ranking_score, 
                length_compound_motif_code, 
                length_protocluster_motif_code,
                highest_ranking_bgc_rank,
                highest_ranking_bgc,
                highest_ranking_bgc_score,
                highest_ranking_bgc_compound_name
            ))

    # sort records on score
    records.sort(key=lambda x: x[2], reverse=True)

    for i, record in enumerate(records):
        (
            compound_id, 
            compound_name, 
            highest_ranking_score, 
            length_compound_motif_code, 
            length_protocluster_motif_code,
            highest_ranking_bgc_rank,
            highest_ranking_bgc,
            highest_ranking_bgc_score,
            highest_ranking_bgc_compound_name
        ) = record
        out_file.write(f"{compound_id},{make_npatlas_url(compound_id)},\"{compound_name}\",{highest_ranking_score},{length_compound_motif_code},{length_protocluster_motif_code},{highest_ranking_bgc_rank},{highest_ranking_bgc},{highest_ranking_bgc_score},{highest_ranking_bgc_compound_name}\n")

    print(f"Number of high scorers: {high_scorers}")
    out_file.close()

    cur.close()
    conn.close()

if __name__ == "__main__":
    main()
