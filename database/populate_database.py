#!/usr/bin/env python3
import argparse 
import atexit
import json
import logging
import os
import typing as ty

import numpy as np
import psycopg2
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from retromol.data.motifs import MOTIFS
from tqdm import tqdm


def setup_logger(logger_level: str) -> logging.Logger:
    """Create a logger with the specified log level, only add a stream handler."""
    # set up logging
    logger = logging.getLogger()
    logger.setLevel(logger_level)

    # add a stream handler
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)

    # ensure the handler is closed properly when exiting
    def close_handler_on_exit():
        for handler in logger.handlers:
            handler.close()
    atexit.register(close_handler_on_exit)

    return logger


def connect_to_database() -> ty.Tuple[psycopg2._ext.connection, psycopg2._ext.cursor]:
    """Connect to the database and return the connection and cursor."""
    # connect to the database
    conn = psycopg2.connect(
        dbname=os.getenv("DB_NAME", "retromol"),
        user=os.getenv("DB_USER", "davidmeijer"),
        password=os.getenv("DB_PASSWORD", "postgres"),
        host=os.getenv("DB_HOST", "localhost"),
        port=os.getenv("DB_PORT", "5432")
    )

    # create cursor and ensure the connection is closed properly when exiting
    cur = conn.cursor()
    def close_connection_on_exit():
        conn.close()
    atexit.register(close_connection_on_exit)

    return conn, cur


def cli() -> argparse.Namespace:
    """Parse command line arguments."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--npatlas", required=False, help="path to the NPAtlas database json file")
    parser.add_argument("--donphan", required=False, help="path to the DONPHAN database csv file")
    parser.add_argument("--retromol", required=False, help="path to RetroMol results folder for NPAtlas compounds")
    parser.add_argument("--log-level", required=False, default="INFO", help="log level for logging", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"])
    return parser.parse_args()


class CompoundRecord:
    def __init__(
        self,
        name: str,
        npa_id: str,
        smiles: str,
        inchikey: str,
    ) -> None:
        if not all(isinstance(arg, str) for arg in (name, npa_id, smiles, inchikey)):
            raise ValueError("all attributes must be strings")
        
        self.name = name
        self.npa_id = npa_id
        self.smiles = smiles
        self.inchikey = inchikey
        self.inchikey_prefix = inchikey.split("-")[0]
    
    def insert_into_database(self, cur: psycopg2._ext.cursor, logger: ty.Optional[logging.Logger] = None) -> None:
        """Insert the compound record into the database."""
        try:
            cur.execute(
                """
                INSERT INTO compounds (name, npa_id, smiles, inchikey, inchikey_prefix)
                VALUES (%s, %s, %s, %s, %s)
                ON CONFLICT (npa_id) DO NOTHING
                """,
                (self.name, self.npa_id, self.smiles, self.inchikey, self.inchikey_prefix)
            )
        except psycopg2.Error as e:
            if logger: logger.error(f"error inserting compound record into database: {e}")


class OrganismRecord:
    def __init__(
            self,
            type: str,
            genus: str,
            species: str,
    ) -> None:
        if not all(isinstance(arg, str) for arg in (type, genus, species)):
            raise ValueError("all attributes must be strings")
        
        self.type = type
        self.genus = genus
        self.species = species

    def insert_into_database(self, cur: psycopg2._ext.cursor, logger: ty.Optional[logging.Logger] = None) -> None:
        """Insert the organism record into the database."""
        try:
            cur.execute(
                """
                INSERT INTO organisms (type, genus, species)
                VALUES (%s, %s, %s)
                ON CONFLICT (genus, species) DO NOTHING
                """,
                (self.type, self.genus, self.species)
            )
        except psycopg2.Error as e:
            if logger: logger.error(f"error inserting organism record into database: {e}")

    def connect_to_compound(self, cur: psycopg2._ext.cursor, compound_npa_id: str, logger: ty.Optional[logging.Logger] = None) -> None:
        """Connect the organism record to a compound record in the database."""
        try:
            cur.execute(
                """
                INSERT INTO compounds_organisms (compound_id, organism_genus, organism_species)
                VALUES (%s, %s, %s)
                ON CONFLICT (compound_id, organism_genus, organism_species) DO NOTHING
                """,
                (compound_npa_id, self.genus, self.species)
            )
        except psycopg2.Error as e:
            if logger: logger.error(f"error connecting organism record to compound record in database: {e}")


class BiosyntheticClassRecord:
    def __init__(
            self,
            name: str,
            level: str,
            classifier: str,
    ) -> None:
        if not all(isinstance(arg, str) for arg in (name, level, classifier)):
            raise ValueError("all attributes must be strings")
        
        self.name = name
        self.level = level
        self.classifier = classifier

    def insert_into_database(self, cur: psycopg2._ext.cursor, logger: ty.Optional[logging.Logger] = None) -> None:
        """Insert the biosynthetic class record into the database."""
        try:
            cur.execute(
                """
                INSERT INTO biosynthetic_classes (name, level, classifier)
                VALUES (%s, %s, %s)
                ON CONFLICT (name, level, classifier) DO NOTHING
                """,
                (self.name, self.level, self.classifier)
            )
        except psycopg2.Error as e:
            if logger: logger.error(f"error inserting biosynthetic class record into database: {e}")

    def connect_to_compound(self, cur: psycopg2._ext.cursor, compound_npa_id: str, logger: ty.Optional[logging.Logger] = None) -> None:
        """Connect the biosynthetic class record to a compound record in the database."""
        try:
            cur.execute(
                """
                INSERT INTO compounds_biosynthetic_classes (compound_id, biosynthetic_class_name, biosynthetic_class_level, biosynthetic_class_classifier)
                VALUES (%s, %s, %s, %s)
                ON CONFLICT (compound_id, biosynthetic_class_name, biosynthetic_class_level, biosynthetic_class_classifier) DO NOTHING
                """,
                (compound_npa_id, self.name, self.level, self.classifier)
            )
        except psycopg2.Error as e:
            if logger: logger.error(f"error connecting biosynthetic class record to compound record in database: {e}")


def connect_bioactivity_to_compound(cur: psycopg2._ext.cursor, compound_npa_id: str, bioactivity_type: str, logger: ty.Optional[logging.Logger] = None) -> None:
    """Connect a bioactivity to a compound record in the database."""
    try:
        cur.execute(
            """
            INSERT INTO compounds_bioactivities (compound_id, bioactivity_type)
            VALUES (%s, %s)
            ON CONFLICT (compound_id, bioactivity_type) DO NOTHING
            """,
            (compound_npa_id, bioactivity_type) 
        )
    except psycopg2.Error as e:
        if logger: logger.error(f"error connecting bioactivity to compound record in database: {e}")


def parse_npatlas(path_to_json: str) -> ty.Generator[ty.Tuple[CompoundRecord, OrganismRecord, ty.List[BiosyntheticClassRecord]], None, None]:
    """Parse the NPAtlas database json file and yield compound records."""
    with open (path_to_json, "r") as f:
        data = json.load(f)
        for record in tqdm(data, desc="parsing NPAtlas database"):

            # parse out compound record
            name = record["original_name"]
            npa_id = record["npaid"]
            smiles = record["smiles"]
            inchikey = record["inchikey"]

            compound_record = CompoundRecord(
                name=name,
                npa_id=npa_id,
                smiles=smiles,
                inchikey=inchikey,
            )
            
            # parse out origin organism record
            origin_organism = record["origin_organism"]
            type_ = origin_organism["type"].lower()
            genus = origin_organism["genus"]
            species = origin_organism["species"]

            organism_record = OrganismRecord(
                type=type_,
                genus=genus,
                species=species,
            )

            # parse out biosynthetic class records
            biosynthetic_class_records = []
            npclassifier_classifications = record.get("npclassifier", {})
            if npclassifier_classifications is None: npclassifier_classifications = {}
            npclassifier_class_results = npclassifier_classifications.get("class_results", [])
            npclassifier_pathway_results = npclassifier_classifications.get("pathway_results", [])
            npclassifier_superclass_results = npclassifier_classifications.get("superclass_results", [])

            for npclassifier_class_result in npclassifier_class_results:
                biosynthetic_class_records.append(
                    BiosyntheticClassRecord(
                        name=npclassifier_class_result,
                        level="class",
                        classifier="npclassifier",
                    )
                )

            for npclassifier_pathway_result in npclassifier_pathway_results:
                biosynthetic_class_records.append(
                    BiosyntheticClassRecord(
                        name=npclassifier_pathway_result,
                        level="pathway",
                        classifier="npclassifier",
                    )
                )

            for npclassifier_superclass_result in npclassifier_superclass_results:
                biosynthetic_class_records.append(
                    BiosyntheticClassRecord(
                        name=npclassifier_superclass_result,
                        level="superclass",
                        classifier="npclassifier",
                    )
                )

            yield compound_record, organism_record, biosynthetic_class_records


def parse_donphan(path_to_csv: str) -> ty.Generator[ty.Tuple[str, ty.List[str]], None, None]:
    """Parse the DONPHAN database csv file and yield bioactivity records."""
    with open(path_to_csv, "r") as f:
        header = f.readline().strip().split(",")[2:]
        for line in tqdm(f, desc="parsing DONPHAN database"):
            # parse out bioactivity record
            _, smiles, *items = line.strip().split(",")
            mol = Chem.MolFromSmiles(smiles)
            inchikey = Chem.inchi.MolToInchiKey(mol)
            inchikey_prefix = inchikey.split("-")[0]

            bioactivity = [len(x) > 0 for x in items]
            items = dict(zip(header, bioactivity))

            # merge keys if they are same but one has '_np' suffix and the other not
            # e.g. 'antibacterial' and 'antibacterial_np'
            for key in list(items.keys()):
                if key.endswith("_np"):
                    new_key = key[:-3]
                    if new_key in items:
                        items[new_key] = items[new_key] or items[key]
            items = {key: value for key, value in items.items() if not key.endswith("_np")}

            # return list of keys where value is True
            bioactivities = [key for key, value in items.items() if value]
        
            yield inchikey_prefix, bioactivities


def smiles_to_fingerprint(smiles: str) -> np.ndarray:
    """Convert a SMILES string to a Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = fp_generator.GetFingerprint(mol)
    fp_arr = np.zeros((2048,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, fp_arr)

    # check if fingerprint only contains 0s and 1s
    assert np.all(np.logical_or(fp_arr == 0, fp_arr == 1))

    return fp_arr

def calc_tanimoto_similarity(arr1: np.ndarray, arr2: np.ndarray) -> float:
    """Calculate the Tanimoto similarity between two fingerprints."""
    # check if arrays are of same size
    assert arr1.shape == arr2.shape

    # check if both arrays only contain 0s and 1s
    assert np.all(np.logical_or(arr1 == 0, arr1 == 1))

    intersection = np.dot(arr1, arr2)
    sum_arr1 = np.sum(arr1)
    sum_arr2 = np.sum(arr2)
    score = intersection / (sum_arr1 + sum_arr2 - intersection)

    # check if score is within [0, 1]
    assert 0 <= score <= 1

    return score


parasect_substrates_as_smiles = {
    "tryptophan": r"C1=CC=C2C(=C1)C(=CN2)C[C@@H](C(=O)O)N",
    "leucine": r"CC(C)C[C@@H](C(=O)O)N",
    "tyrosine": r"C1=CC(=CC=C1C[C@@H](C(=O)O)N)O",
    "serine": r"C([C@@H](C(=O)O)N)O",
    "alanine": r"C[C@@H](C(=O)O)N",
    "valine": r"CC(C)[C@@H](C(=O)O)N",
    "arginine": r"C(C[C@@H](C(=O)O)N)CN=C(N)N",
    "glutamine": r"C(CC(=O)N)[C@@H](C(=O)O)N",
    "isoleucine": r"CC[C@H](C)[C@@H](C(=O)O)N",
    "aspartic acid": r"C([C@@H](C(=O)O)N)C(=O)O",
    "cysteine": r"C([C@@H](C(=O)O)N)S",
    "proline": r"C1C[C@H](NC1)C(=O)O",
    "phenylalanine": r"C1=CC=C(C=C1)C[C@@H](C(=O)O)N",
    "asparagine": r"C([C@@H](C(=O)O)N)C(=O)N",
    "4-hydroxyphenylglycine": r"C1=CC(=CC=C1[C@@H](C(=O)O)N)O",
    "2,4-diaminobutyric acid": r"C(CN)[C@@H](C(=O)O)N",
    "beta-hydroxytyrosine": r"C1=CC(=CC=C1[C@H]([C@@H](C(=O)O)N)O)O",
    "glutamic acid": r"C(CC(=O)O)[C@@H](C(=O)O)N",
    "glycine": r"NCC(=O)O",
    "3,5-dihydroxyphenylglycine": r"N[C@H](C(=O)O)c1cc(O)cc(O)c1",
    "pipecolic acid": r"C1CCN[C@@H](C1)C(=O)O",
    "threonine": r"C[C@H]([C@@H](C(=O)O)N)O",
    "N5-formyl-N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CN(C=O)O",
    "N5-hydroxyornithine": r"C(C[C@@H](C(=O)O)N)CNO",
    "ornithine": r"C(C[C@@H](C(=O)O)N)CN",
    "beta-alanine": r"NCCC(=O)O",
    "histidine": r"C1=C(NC=N1)C[C@@H](C(=O)O)N",
    "D-alanine": r"C[C@H](C(=O)O)N",
    "2-aminoadipic acid": r"C(C[C@@H](C(=O)O)N)CC(=O)O",
    "lysine": r"C(CCN)C[C@@H](C(=O)O)N",
    "2-aminoisobutyric acid": r"O=C(O)C(N)(C)C",
    "2,3-dihydroxybenzoic acid": r"C1=CC(=C(C(=C1)O)O)C(=O)O",
    "salicylic acid": r"C1=CC=C(C(=C1)C(=O)O)O",
    "anthranillic acid": r"C1=CC=C(C(=C1)C(=O)O)N",
}


polyketide_motif_identities = [
    "A", "B", "C", "D",
    "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11",
    "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11",
    "C1", "C2", "C4",
    "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
]


# convert parasect substrates to fingerprints
substrate_to_fingerprint = {
    name: smiles_to_fingerprint(smiles) 
    for name, smiles 
    in tqdm(parasect_substrates_as_smiles.items(), desc="converting parasect substrates to fingerprints...", total=len(parasect_substrates_as_smiles), leave=True)
}

# convert retromol motifs to fingerprints
motif_to_fingerprint = {
    motif["name"]: smiles_to_fingerprint(motif["smiles"])
    for motif in tqdm(MOTIFS, desc="converting retromol motifs to fingerprints...", total=len(MOTIFS), leave=True)
}


def format_polyketide_tautomer_identity(name: str) -> str:
    if "-tautomer" in name:
        return name.replace("-tautomer", "")
    return name


polyketide_smarts = Chem.MolFromSmarts("[OH]SC~CC(=O)[OH]")
alpha_amino_acid_smarts = Chem.MolFromSmarts("[NH2]CC(=O)[OH]")
alpha_amino_acid_smart_proline_like = Chem.MolFromSmarts("[NH1]CC(=O)[OH]")
beta_amino_acid_smarts = Chem.MolFromSmarts("[NH2]CCC(=O)[OH]")


def match_other_motif(smiles: str) -> ty.Optional[str]:
    mol = Chem.MolFromSmiles(smiles)

    if mol.HasSubstructMatch(polyketide_smarts):
        return "polyketide"
    if mol.HasSubstructMatch(alpha_amino_acid_smarts):
        return "amino acid"
    if mol.HasSubstructMatch(alpha_amino_acid_smart_proline_like):
        return "amino acid"
    if mol.HasSubstructMatch(beta_amino_acid_smarts):
        return "amino acid"
    
    return None


def parse_motif(name: ty.Optional[str], smiles: str) -> str:
    # unidentified motifs are named 'other'
    if name is None or name == "other":
        return "other"

    # check if start or end
    if name == "start" or name == "end":
        return name

    # check for '-tautomer' suffix in polyketide motif identities
    name = format_polyketide_tautomer_identity(name)

    # parse name
    if name in polyketide_motif_identities:
        return name

    # check if it is substrate known by parasect, check by comparing fingerprints
    motif_fp = motif_to_fingerprint[name]
    for substrate_name, substrate_fp in substrate_to_fingerprint.items():
        score = calc_tanimoto_similarity(motif_fp, substrate_fp)
        if score == 1.0:
            return substrate_name
        
    # check if we can identify a backbone in the motif
    if name := match_other_motif(smiles):
        return name
    
    # unidentified motifs are named 'other'
    return "other"


def parse_retromol_results(results_folder: str) -> ty.Generator[ty.Tuple[str, ty.List[str]], None, None]:
    # parse over subfolders, select for put.json file in each subfolder
    for subdir in tqdm(os.listdir(results_folder), desc="parsing retromol results", leave=True):
        data_path = os.path.join(results_folder, subdir, "out.json")

        # check if file exists, if not the item failed to parse
        if not os.path.exists(data_path):
            continue

        with open(data_path, "r") as f:
            data = json.load(f)

        # parse data
        npa_id = data["name"]
        encoding_to_motif_codes = data["encoding_to_motif_codes"]

        # parse out sequences
        for _, motif_codes in encoding_to_motif_codes.items():
            for motif_code in motif_codes.values():
                primary_sequence = ["start"]
                for motif in motif_code[::-1]:
                    motif_name = motif["identity"]
                    motif_smiles = motif["smiles"]
                    parsed_motif_string = parse_motif(motif_name, motif_smiles)
                    primary_sequence.append(parsed_motif_string)
                primary_sequence.append("end")

                yield npa_id, primary_sequence


def main() -> None:
    """Main function."""
    # suppress RDKit warnings
    RDLogger.DisableLog("rdApp.*")

    # parse command line arguments
    args = cli()
    logger_level = args.log_level
    npatlas_path = args.npatlas
    donphan_path = args.donphan
    retromol_path = args.retromol

    # set up logging
    logger = setup_logger(logger_level)
    logger.info("populating the database...")

    # create database connection and cursor
    conn, cur = connect_to_database()

    # populate the database with NPAtlas compounds
    if npatlas_path is not None:
        for compound_record, organism_record, biosynthetic_class_records in parse_npatlas(npatlas_path):
            compound_record.insert_into_database(cur, logger)
            organism_record.insert_into_database(cur, logger)   
            organism_record.connect_to_compound(cur, compound_record.npa_id, logger)
            for biosynthetic_class_record in biosynthetic_class_records:
                biosynthetic_class_record.insert_into_database(cur, logger)
                biosynthetic_class_record.connect_to_compound(cur, compound_record.npa_id, logger)
        conn.commit()

    # count the number of compound records in the database
    cur.execute("SELECT COUNT(*) FROM compounds")
    count = cur.fetchone()[0]
    logger.info(f"number of compound records in the database: {count}")

    # count the number of organism records in the database
    cur.execute("SELECT COUNT(*) FROM organisms")
    count = cur.fetchone()[0]
    logger.info(f"number of organism records in the database: {count}")

    # count the number of biosynthetic class records in the database
    cur.execute("SELECT COUNT(*) FROM biosynthetic_classes")
    count = cur.fetchone()[0]
    logger.info(f"number of biosynthetic_class records in the database: {count}")

    # count the number of bioactivity records in the database
    cur.execute("SELECT COUNT(*) FROM bioactivities")
    count = cur.fetchone()[0]
    logger.info(f"number of bioactivity records in the database: {count}")

    # count the number of motifs records in the database
    cur.execute("SELECT COUNT(*) FROM motifs")
    count = cur.fetchone()[0]
    logger.info(f"number of motifs records in the database: {count}")

    # populate the database with DONPHAN bioactivities
    if donphan_path is not None:

        # retrieve all compound inchikey prefixes and npa_ids from the database
        compound_inchikey_prefix_to_npa_ids = {}
        cur.execute("SELECT npa_id, inchikey_prefix FROM compounds")
        for npa_id, inchikey_prefix in cur.fetchall():
            compound_inchikey_prefix_to_npa_ids.setdefault(inchikey_prefix, []).append(npa_id)

        # parse the DONPHAN database
        for inchikey_prefix, bioactivities in parse_donphan(donphan_path):
            npa_ids = compound_inchikey_prefix_to_npa_ids.get(inchikey_prefix, [])
            for npa_id in npa_ids:
                for bioactivity in bioactivities:
                    connect_bioactivity_to_compound(cur, npa_id, bioactivity, logger)
        conn.commit()

    # populate the database with RetroMol results
    if retromol_path is not None:
        retromol_version = "1.0.0-dev"

        for npa_id, primary_sequence in parse_retromol_results(retromol_path):
            # insert into primary_sequences
            cur.execute(
                """
                INSERT INTO primary_sequences (retromol_version)
                VALUES (%s)
                RETURNING id;
                """,
                (retromol_version,)
            )
            primary_sequence_id = cur.fetchone()[0]

            # link the sequence to the compound
            cur.execute(
                """
                INSERT INTO compounds_primary_sequences (compound_id, primary_sequence_id) 
                VALUES (%s, %s);
                """,
                (npa_id, primary_sequence_id)
            )

            # add motifs for the sequence
            for position, motif in enumerate(primary_sequence):
                cur.execute(
                    """
                    INSERT INTO primary_sequences_motifs (primary_sequence_id, position, motif_id)
                    VALUES (%s, %s, %s);
                    """,
                    (primary_sequence_id, position, motif)
                )
        conn.commit()

    # count the number of primary_sequence records in the database
    cur.execute("SELECT COUNT(*) FROM primary_sequences")
    count = cur.fetchone()[0]
    logger.info(f"number of primary_sequence records in the database: {count}")

    # close cursor and connection
    cur.close()
    conn.close()

    # TODO: parse and link antismash database records
    # TODO: parse and link mibig database records to antismash records (download, parse, and put in next to the antismash database records)

    logger.info("done")


if __name__ == "__main__":
    main()
