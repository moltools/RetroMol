import psycopg2
from rdkit import Chem, RDLogger
from tqdm import tqdm


def parse_donphan(path_donphan):
    """Parse the DONPHAN database csv file and yield bioactivity records."""
    with open(path_donphan, "r") as f:
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


def connect_compound_to_bioactivity(cur, compound_id, bioactivity):
    """Connect a compound to a bioactivity."""
    try:
        cur.execute(
            """
            INSERT INTO compounds_bioactivities (compound_id, bioactivity_type)
            VALUES (%s, %s)
            ON CONFLICT (compound_id, bioactivity_type) DO NOTHING
            """,
            (compound_id, bioactivity)
        )
    except psycopg2.Error as e:
        raise ValueError(f"error connecting compound record to bioactivity record in database: {e}")


def add_donphan(cur, path_donphan):
    """Add bioactivity records from the DONPHAN database to the database."""
    RDLogger.DisableLog("rdApp.*")

    # retrieve all compound inchikey prefixes and npa_ids from the database
    compound_inchikey_prefix_to_compound_id = {}
    cur.execute("SELECT compound_id, inchikey_prefix FROM compounds")
    for compound_id, inchikey_prefix in cur.fetchall():
        compound_inchikey_prefix_to_compound_id.setdefault(inchikey_prefix, []).append(compound_id)

    # parse DONPHAN database file
    for inchikey_prefix, bioactivities in parse_donphan(path_donphan):
        compound_ids = compound_inchikey_prefix_to_compound_id.get(inchikey_prefix, [])
        for compound_id in compound_ids:
            for bioactivity in bioactivities:
                connect_compound_to_bioactivity(cur, compound_id, bioactivity)
