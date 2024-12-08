#!/usr/bin/env python3
import argparse
import json
import os
from typing import Dict, List

from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem
from tqdm import tqdm


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mibig", required=True, type=str, help="path to mibig jsons folder")
    parser.add_argument("--dataset", required=True, type=str, help="path to dataset csv file")
    parser.add_argument("--output", required=True, type=str, help="path to output folder")
    return parser.parse_args()


def parse_datset(path_to_file: str) -> Dict[str, str]:
    """Parses combpounds from dataset file in dictionary as {compound_id: inchikey_connectivity, ...}."""
    with open(path_to_file, "r") as f:
        f.readline() # skip header
        lines = f.readlines()
        dataset = {}
        for line in tqdm(lines):
            compound_id, smiles, *_ = line.strip().split(",")
            mol = Chem.MolFromSmiles(smiles)
            inchikey_connectivity = Chem.MolToInchiKey(mol).split("-")[0]
            dataset[compound_id] = inchikey_connectivity
    return dataset


class MibigEntry:
    def __init__(self, accession: str, compound_connectivities: List[str]) -> None:
        self.accession = accession
        self.compound_connectivities = compound_connectivities


def parse_mibig(path_to_folder: str) -> List[MibigEntry]:
    """Parses mibig jsons from folder."""
    records = []
    # loop over all jsons in folder
    for json_file in tqdm(os.listdir(path_to_folder)):
        if not json_file.endswith(".json"):
            continue
        json_file_path = os.path.join(path_to_folder, json_file)
        with open(os.path.join(path_to_folder, json_file), "r") as f:
            mibig_entry = json.load(f)
            accession = mibig_entry["accession"]
            compound_inchikey_connectivities = []
            for compound_entry in mibig_entry.get("compounds", []):
                compound_smiles = compound_entry.get("structure", None)
                if compound_smiles is not None:
                    compound_mol = Chem.MolFromSmiles(compound_smiles)
                    compound_inchikey_connectivity = Chem.MolToInchiKey(compound_mol).split("-")[0]
                    compound_inchikey_connectivities.append(compound_inchikey_connectivity)
            if len(compound_inchikey_connectivities) > 0:
                records.append(MibigEntry(accession, compound_inchikey_connectivities))
    return records


def main() -> None:
    RDLogger.DisableLog("rdApp.*")
    args = cli()    
    included_compounds = parse_datset(args.dataset)
    print("Number of included compounds:", len(included_compounds))
    mibig_entries = parse_mibig(args.mibig)
    print("Number of mibig entries with compounds:", len(mibig_entries))
    mibig_accessions_with_compounds_in_dataset = []
    for mibig_entry in mibig_entries:
        if any([compound_connectivity in included_compounds.values() for compound_connectivity in mibig_entry.compound_connectivities]):
            mibig_accessions_with_compounds_in_dataset.append(mibig_entry.accession)
    print("Number of mibig entries with compounds in dataset:", len(mibig_accessions_with_compounds_in_dataset))
    mibig_accessions_with_compounds_in_dataset.sort()
    with open(os.path.join(args.output, "mibig_accessions_with_compounds_in_dataset.txt"), "w") as f:
        for accession in mibig_accessions_with_compounds_in_dataset:
            f.write(accession + "\n")


if __name__ == "__main__":
    main()