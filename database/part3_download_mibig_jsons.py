#!/usr/bin/env python3
import argparse 
import json
import os
import requests
import time
import typing as ty

from rdkit import Chem, RDLogger
from tqdm import tqdm

from _utils import connect_to_database, setup_logger
from part2_select_npatlas_compounds import retrieve_compounds


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mibig", type=str, required=True, help="path to mibig jsons folder")
    parser.add_argument("--outdir", type=str, required=True, help="path to output dir")
    parser.add_argument("--log-level", required=False, default="INFO", help="log level for logging", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args()


class MibigEntry:
    def __init__(self, accession, version, compound_connectivities) -> None:
        self._accession = accession
        self._version = version
        self.compound_connectivities = compound_connectivities

    @property
    def accession(self):
        return f"{self._accession}.{self._version}"


def parse_mibig(path_to_folder: str, logger) -> ty.List[MibigEntry]:
    """Parses mibig jsons from folder."""
    records = []
    for json_file in tqdm(os.listdir(path_to_folder)):

        if not json_file.endswith(".json"):
            continue

        with open(os.path.join(path_to_folder, json_file), "r") as f:
            mibig_entry = json.load(f)
            accession = mibig_entry["accession"]
            accession_version = mibig_entry["version"]
            
            compound_inchikey_connectivities = []
            for compound_entry in mibig_entry.get("compounds", []):
                compound_smiles = compound_entry.get("structure", None)
                if compound_smiles is not None:
                    compound_mol = Chem.MolFromSmiles(compound_smiles)
                    compound_inchikey_connectivity = Chem.MolToInchiKey(compound_mol).split("-")[0]
                    compound_inchikey_connectivities.append(compound_inchikey_connectivity)

            if len(compound_inchikey_connectivities) > 0:
                records.append(MibigEntry(accession, accession_version, compound_inchikey_connectivities))

    logger.info(f"Found {len(records)} MIBiG entries")
    return records


def download_single_bgc(bgc, save_path, logger):
    """Download a single BGC from MIBiG."""
    base_url = "https://mibig.secondarymetabolites.org/repository/{}/generated/{}.json"

    try:
        # Construct the URL
        url = base_url.format(bgc, bgc.split('.')[0])
        response = requests.get(url)
        
        if response.status_code == 200:
            with open(save_path, 'w', encoding='utf-8') as f:
                f.write(response.text)
            logger.info(f"Downloaded: {save_path}")
        else:
            logger.error(f"Failed to download {bgc} from {url}. HTTP Status Code: {response.status_code}")
            logger.error(f"BGC with accession {bgc} does not exist in MIBiG database and is probably retired.")
    except Exception as e:
        logger.error(f"An error occurred for {bgc}: {e}")
    
    # Wait for 1 second to avoid overloading the server
    time.sleep(1)


def main() -> None:
    # disable rdkit logging
    RDLogger.DisableLog("rdApp.*")

    # parse arguments
    args = cli()
    log_level = args.log_level
    outdir = args.outdir
    mibig_folder = args.mibig

    # setup logger
    logger = setup_logger(log_level)
    logger.info("Selecting NPAtlas compounds")

    # parse mibig entries
    mibig_records = parse_mibig(mibig_folder, logger)
    accession_to_record = {r.accession: r for r in mibig_records}

    # retrieve npatlas compounds
    conn, cur = connect_to_database()
    npatlas_compounds = retrieve_compounds(cur, logger, return_inchikey=True)
    conn.close()
    cur.close()
    included_compounds = {v.split("-")[0]: k for k, v in npatlas_compounds}  # {inchikey_suffix: compound_id, ...}

    # find ouw which mibig entries have compounds in npatlas dataset
    mibig_accessions_with_compounds_in_dataset = {}
    for mibig_record in mibig_records:
        for compound_connectivity in mibig_record.compound_connectivities:
            if compound_connectivity in included_compounds:
                if mibig_record.accession not in mibig_accessions_with_compounds_in_dataset:
                    mibig_accessions_with_compounds_in_dataset[mibig_record.accession] = [included_compounds[compound_connectivity]]
                else:
                    mibig_accessions_with_compounds_in_dataset[mibig_record.accession].append(included_compounds[compound_connectivity])

    print("Number of mibig entries with compounds in dataset:", len(mibig_accessions_with_compounds_in_dataset))
    print("Unique mibig accessions with compounds in dataset:", len(set(mibig_accessions_with_compounds_in_dataset)))
    
    # check files int outfolder and list the accession that are still missing
    downloaded_files = os.listdir(outdir)
    missing_accessions = []
    for accession in mibig_accessions_with_compounds_in_dataset.keys():
        if accession.split(".")[0] + ".json" not in downloaded_files:
            missing_accessions.append(accession)
    print("Number of missing accessions:", len(missing_accessions))

    if len(missing_accessions) > 0:
        accessions_to_download = missing_accessions
    else:
        accessions_to_download = list(mibig_accessions_with_compounds_in_dataset.keys())

    # download mibig jsons
    for accession in accessions_to_download:
        save_path = os.path.join(outdir, accession.split(".")[0] + ".json")
        download_single_bgc(accession, save_path, logger)

    logger.info("Finished downloading MIBiG jsons")


if __name__ == "__main__":
    main()
