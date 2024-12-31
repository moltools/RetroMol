import json
import os
import numpy as np
from tqdm import tqdm
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from retromol.data.motifs import MOTIFS


RETROMOL_VERSION = "1.0.0-dev"


def smiles_to_fingerprint(smiles):
    """Convert a SMILES string to a Morgan fingerprint."""
    mol = Chem.MolFromSmiles(smiles, sanitize=True)
    fp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fp = fp_generator.GetFingerprint(mol)
    fp_arr = np.zeros((2048,), dtype=np.int8)
    DataStructs.ConvertToNumpyArray(fp, fp_arr)
    assert np.all(np.logical_or(fp_arr == 0, fp_arr == 1))
    return fp_arr

def calc_tanimoto_similarity(arr1, arr2):
    """Calculate the Tanimoto similarity between two fingerprints."""
    assert arr1.shape == arr2.shape
    assert np.all(np.logical_or(arr1 == 0, arr1 == 1))
    intersection = np.dot(arr1, arr2)
    sum_arr1 = np.sum(arr1)
    sum_arr2 = np.sum(arr2)
    score = intersection / (sum_arr1 + sum_arr2 - intersection)
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
    "R-beta-hydroxytyrosine": r"C1=CC(=CC=C1[C@H]([C@@H](C(=O)O)N)O)O",
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
    "anthranilic acid": r"C1=CC=C(C(=C1)C(=O)O)N",
}


polyketide_motif_identities = [
    "A", "B", "C", "D",
    "A1", "A2", "A3", "A4", "A5", "A6", "A7", "A8", "A9", "A10", "A11",
    "B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10", "B11",
    "C1", "C2", "C4",
    "D1", "D2", "D3", "D4", "D5", "D6", "D7", "D8", "D9", "D10", "D11", "D12",
]


# convert parasect substrates to fingerprints
substrate_to_fingerprint = {name: smiles_to_fingerprint(smiles) for name, smiles in parasect_substrates_as_smiles.items()}
# convert retromol motifs to fingerprints
motif_to_fingerprint = {motif["name"]: smiles_to_fingerprint(motif["smiles"]) for motif in MOTIFS}

polyketide_smarts = Chem.MolFromSmarts("[OH]SC~CC(=O)[OH]")
alpha_amino_acid_smarts = Chem.MolFromSmarts("[NH2]CC(=O)[OH]")
alpha_amino_acid_smart_proline_like = Chem.MolFromSmarts("[NH1]CC(=O)[OH]")
beta_amino_acid_smarts = Chem.MolFromSmarts("[NH2]CCC(=O)[OH]")


def match_other_motif(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol.HasSubstructMatch(polyketide_smarts): return "polyketide"
    if mol.HasSubstructMatch(alpha_amino_acid_smarts): return "amino acid"
    if mol.HasSubstructMatch(alpha_amino_acid_smart_proline_like): return "amino acid"
    if mol.HasSubstructMatch(beta_amino_acid_smarts): return "amino acid"
    return None


def parse_motif(name, smiles):
    """Translate primary sequence motif to a more readable format."""
    # unidentified motifs are named 'other'
    if name is None or name == "other": return "other"

    # check if start or end
    if name == "start" or name == "end": return name

    # check for '-tautomer' suffix in polyketide motif identities
    if "-tautomer" in name: name = name.replace("-tautomer", "")

    # parse name
    if name in polyketide_motif_identities: return name

    # check if it is substrate known by parasect, check by comparing fingerprints
    motif_fp = motif_to_fingerprint[name]
    for substrate_name, substrate_fp in substrate_to_fingerprint.items():
        score = calc_tanimoto_similarity(motif_fp, substrate_fp)
        if score == 1.0: return substrate_name
        
    # check if we can identify a backbone in the motif
    if name := match_other_motif(smiles): return name
    
    # unidentified motifs are named 'other'
    return "other"


def parse_retromol_results(path_retromol_results):
    """Parse primary sequences from RetroMol results folder."""
    RDLogger.DisableLog("rdApp.*")

    # parse over subfolders, select for put.json file in each subfolder
    for subdir in tqdm(os.listdir(path_retromol_results), desc="parsing retromol results", leave=True):
        data_path = os.path.join(path_retromol_results, subdir, "out.json")

        # check if file exists, if not the item failed to parse
        if not os.path.exists(data_path): continue

        with open(data_path, "r") as f: data = json.load(f)

        # parse data
        compound_id = data["name"]
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
                yield compound_id, primary_sequence


def connect_compound_to_primary_sequence(cur, compound_source, compound_id, primary_sequence):
    """Add primary sequence to database and connect to compound it belongs to."""
    # add new primary sequence
    cur.execute(
        """
        INSERT INTO primary_sequences (retromol_version)
        VALUES (%s)
        RETURNING id;
        """,
        (RETROMOL_VERSION,)
    )
    primary_sequence_id = cur.fetchone()[0]

    # link primary sequence to compound it originates from
    cur.execute(
        """
        INSERT INTO compounds_primary_sequences (compound_id, primary_sequence_id)
        VALUES (%s, %s);
        """,
        (compound_id, primary_sequence_id)
    )

    # add primary sequence motifs
    for position, motif in enumerate(primary_sequence):
        cur.execute(
            """
            INSERT INTO primary_sequences_motifs (primary_sequence_id, position, motif_id)
            VALUES (%s, %s, %s);
            """,
            (primary_sequence_id, position, motif)
        )


def add_retromol_results(cur, path_retromol_results, compound_source="npatlas"):
    """Add Retromol results to the database."""
    for compound_id, primary_sequence in parse_retromol_results(path_retromol_results):
        connect_compound_to_primary_sequence(cur, compound_source, compound_id, primary_sequence)
