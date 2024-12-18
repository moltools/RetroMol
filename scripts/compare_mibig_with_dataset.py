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
    def __init__(self, accession: str, compound_connectivities: List[str], primary_sequence: List[str]):
        self.accession = accession
        self.compound_connectivities = compound_connectivities
        self.primary_sequence = primary_sequence


def parse_mibig(path_to_folder: str) -> List[MibigEntry]:
    """Parses mibig jsons from folder."""
    records = []
    # loop over all jsons in folder
    for json_file in tqdm(os.listdir(path_to_folder)):
        if not json_file.endswith(".json"):
            continue
        with open(os.path.join(path_to_folder, json_file), "r") as f:
            mibig_entry = json.load(f)
            accession = mibig_entry["accession"]
            biosynthesis = mibig_entry["biosynthesis"]
            primary_sequence = []
            if modules := biosynthesis.get("modules", None):
                for module in modules:
                    module_type = module["type"]
                    if domains := module.get("modification_domains", None): module_modifications = [domain["type"] for domain in domains]
                    else: module_modifications = []
                    if module_type == "nrps-type1":
                        primary_sequence.append("amino acid")
                    elif module_type == "pks-modular":
                        if accession == "BGC0000055": print(module_modifications)
                        if "ketoreductase" in module_modifications and "dehydratase" in module_modifications and "enoylreductase" in module_modifications:
                            primary_sequence.append("D")
                        elif "ketoreductase" in module_modifications and "dehydratase" in module_modifications:
                            primary_sequence.append("C")
                        elif "ketoreductase" in module_modifications:
                            primary_sequence.append("B")
                        else:
                            primary_sequence.append("A")
                    elif module_type == "pks-trans-at-starter" or module_type == "pks-modular-starter":
                        primary_sequence.append("other")
            compound_inchikey_connectivities = []
            for compound_entry in mibig_entry.get("compounds", []):
                compound_smiles = compound_entry.get("structure", None)
                if compound_smiles is not None:
                    compound_mol = Chem.MolFromSmiles(compound_smiles)
                    compound_inchikey_connectivity = Chem.MolToInchiKey(compound_mol).split("-")[0]
                    compound_inchikey_connectivities.append(compound_inchikey_connectivity)
            if len(compound_inchikey_connectivities) > 0:
                records.append(MibigEntry(accession, compound_inchikey_connectivities, primary_sequence))
    return records


def main() -> None:
    RDLogger.DisableLog("rdApp.*")
    args = cli()    
    included_compounds = parse_datset(args.dataset)  # {compound_id: inchikey_connectivity, ...}
    included_compounds = {v: k for k, v in included_compounds.items()}  # {inchikey_connectivity: compound_id, ...}
    print("Number of included compounds:", len(included_compounds))

    mibig_entries = parse_mibig(args.mibig)
    accession_to_record = {entry.accession: entry for entry in mibig_entries}
    print("Number of mibig entries with compounds:", len(mibig_entries))

    mibig_accessions_with_compounds_in_dataset = {}
    for mibig_entry in mibig_entries:
        for compound_connectivity in mibig_entry.compound_connectivities:
            if compound_connectivity in included_compounds:
                if mibig_entry.accession not in mibig_accessions_with_compounds_in_dataset:
                    mibig_accessions_with_compounds_in_dataset[mibig_entry.accession] = [included_compounds[compound_connectivity]]
                else:
                    mibig_accessions_with_compounds_in_dataset[mibig_entry.accession].append(included_compounds[compound_connectivity])

    print("Number of mibig entries with compounds in dataset:", len(mibig_accessions_with_compounds_in_dataset))
    print("Unique mibig accessions with compounds in dataset:", len(set(mibig_accessions_with_compounds_in_dataset)))

    mibig_accession_keys = sorted(list(mibig_accessions_with_compounds_in_dataset.keys()))
    with open(os.path.join(args.output, "mibig_accessions_with_compounds_in_dataset.csv"), "w") as f:
        f.write("item,mibig,primary_sequence,npaids\n")
        for i, accession in enumerate(mibig_accession_keys):
            primary_sequence = "|".join(accession_to_record[accession].primary_sequence)
            f.write(f"{i+1},{accession},{primary_sequence},{'|'.join(sorted(mibig_accessions_with_compounds_in_dataset[accession]))}\n")


if __name__ == "__main__":
    main()