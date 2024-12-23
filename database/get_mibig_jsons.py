#!/usr/bin/env python3
import argparse 
import json
import os
import requests
import time
from tqdm import tqdm


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mibig-database", type=str, required=True, help="path to mibig database folder")
    parser.add_argument("--outdir", type=str, required=True, help="path to output dir")
    return parser.parse_args()


def parse_mibig(path_mibig_database):
    """Parse mibig accession strings from mibig database."""
    for json_file in tqdm(os.listdir(path_mibig_database)):
        if not json_file.endswith(".json"): continue
        with open(os.path.join(path_mibig_database, json_file), "r") as f:
            mibig_entry = json.load(f)
            accession = mibig_entry["accession"]
            accession_version = mibig_entry["version"]
            yield accession, accession_version


def download_mibig_json(bgc, save_path):
    """Download a single BGC from MIBiG."""
    base_url = "https://mibig.secondarymetabolites.org/repository/{}/generated/{}.json"
    try:
        url = base_url.format(bgc, bgc.split('.')[0])
        response = requests.get(url)
        if response.status_code == 200:
            with open(save_path, 'w', encoding='utf-8') as f: f.write(response.text)
        else: print(f"Failed to download {bgc} from {url}. HTTP Status Code: {response.status_code}. Entry might be retired.")
    except Exception as e: print(f"Failed to download {bgc} from {url}. Error: {e}")
    time.sleep(0.1)


def main() -> None:
    args = cli()
    path_mibig_database = args.mibig_database
    path_outdir = args.outdir

    # get mibig accesssions
    mibig_accessions = [f"{accession}.{version}" for accession, version in parse_mibig(path_mibig_database)]
    print(f"number of MIBiG entries: {len(mibig_accessions)}")

    # check which mibig accessions are already downloaded and inside path_outdir
    downloaded_files = os.listdir(path_outdir)
    missing_accessions = [accession for accession in mibig_accessions if accession.split(".")[0] + ".json" not in downloaded_files]
    print(f"number of missing MIBiG entries up for downloading: {len(missing_accessions)}")

    # download mibig jsons
    for accession in tqdm(missing_accessions):
        save_path = os.path.join(path_outdir, accession.split(".")[0] + ".json")
        download_mibig_json(accession, save_path)


if __name__ == "__main__":
    main()
