import json
import psycopg2
from tqdm import tqdm


class CompoundRecord:
    def __init__(self, name, compound_id, smiles, inchikey):
        if not all(isinstance(arg, str) for arg in (name, compound_id, smiles, inchikey)):
            raise ValueError("all attributes must be strings")
        
        self.name = name
        self.compound_id = compound_id
        self.smiles = smiles
        self.inchikey = inchikey
        self.inchikey_prefix = inchikey.split("-")[0]
    
    def insert_into_database(self, cur):
        """Insert the compound record into the database."""
        try:
            cur.execute(
                """
                INSERT INTO compounds (compound_id, source, name, smiles, inchikey, inchikey_prefix)
                VALUES (%s, %s, %s, %s, %s, %s)
                ON CONFLICT (compound_id) DO NOTHING
                """,
                (self.compound_id, "npatlas", self.name, self.smiles, self.inchikey, self.inchikey_prefix)
            )
        except psycopg2.Error as e:
            raise ValueError(f"error inserting compound record into database: {e}")


class OrganismRecord:
    def __init__(self, type_, genus, species):
        if not all(isinstance(arg, str) for arg in (type_, genus, species)):
            raise ValueError("all attributes must be strings")
        
        assert type_ in ["bacterium", "fungus", "plant", "unknown"], f"invalid organism type: {_type}"
        
        # if genus contains 'unknown', set both genus and species to None
        if genus is not None and species is not None:
            if "unknown" in genus.lower():
                genus = "unknown"
                species = "unknown"

        # if only species contains 'unknown' or 'sp.', set species to None
        if species is not None:
            if "unknown" in species.lower() or "sp." in species.lower():
                species = "unknown"
        
        self.type_ = type_
        self.genus = genus
        self.species = species

    def insert_into_database(self, cur):
        """Insert the organism record into the database."""
        try:
            cur.execute(
                """
                INSERT INTO organisms (type, genus, species)
                VALUES (%s, %s, %s)
                ON CONFLICT (type, genus, species)
                DO UPDATE SET type = organisms.type
                RETURNING id
                """,
                (self.type_, self.genus, self.species)
            )
            organism_id = cur.fetchone()[0]
            return organism_id
        except psycopg2.Error as e:
            raise ValueError(f"error inserting organism record into database: {e}")


class BiosyntheticClassRecord:
    def __init__(self, name, level, classifier):
        if not all(isinstance(arg, str) for arg in (name, level, classifier)):
            raise ValueError("all attributes must be strings")
        
        self.name = name
        self.level = level
        self.classifier = classifier

    def insert_into_database(self, cur) -> None:
        """Insert the biosynthetic class record into the database."""
        try:
            cur.execute(
                """
                INSERT INTO biosynthetic_classes (name, level, classifier)
                VALUES (%s, %s, %s)
                ON CONFLICT (name, level, classifier)
                DO UPDATE SET name = biosynthetic_classes.name
                RETURNING id
                """,
                (self.name, self.level, self.classifier)
            )
            biosynthetic_class_id = cur.fetchone()[0]
            return biosynthetic_class_id
        except psycopg2.Error as e:
            raise ValueError(f"error inserting biosynthetic class record into database: {e}")


def connect_compound_to_organism(cur, compound_id, organism_id) -> None:
    """Connect the organism record to a compound record in the database."""
    try:
        cur.execute(
            """
            INSERT INTO compounds_organisms (compound_id, organism_id)
            VALUES (%s, %s)
            ON CONFLICT (compound_id, organism_id) DO NOTHING
            """,
            (compound_id, organism_id)
        )
    except psycopg2.Error as e:
        raise ValueError(f"error connecting organism record to compound record in database: {e}")


def connect_compound_to_biosynthetic_class(cur, compound_id, biosynthetic_class_id) -> None:
    """Connect the biosynthetic class record to a compound record in the database."""
    try:
        cur.execute(
            """
            INSERT INTO compounds_biosynthetic_classes (compound_id, biosynthetic_class_id)
            VALUES (%s, %s)
            ON CONFLICT (compound_id, biosynthetic_class_id) DO NOTHING
            """,
            (compound_id, biosynthetic_class_id)
        )
    except psycopg2.Error as e:
        raise ValueError(f"error connecting biosynthetic class record to compound record in database: {e}")


def parse_npatlas(path_npatlas):
    """Parse the NPAtlas database json file and yield compound records."""
    with open (path_npatlas, "r") as f:
        data = json.load(f)
        for record in tqdm(data, desc="parsing NPAtlas database"):

            # parse out compound record
            name = record["original_name"]
            npaid = record["npaid"]
            smiles = record["smiles"]
            inchikey = record["inchikey"]

            compound_record = CompoundRecord(
                name=name,
                compound_id=npaid,
                smiles=smiles,
                inchikey=inchikey,
            )
            
            # parse out origin organism record
            origin_organism = record["origin_organism"]
            type_ = origin_organism["type"].lower()
            genus = origin_organism["genus"]
            species = origin_organism["species"]

            organism_record = OrganismRecord(
                type_=type_,
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


def add_npatlas(cur, path_npatlas):
    """Add NPAtlas database records to the database."""
    for compound_record, organism_record, biosynthetic_class_records in parse_npatlas(path_npatlas):
        # insert compound record into database
        compound_record.insert_into_database(cur)

        # insert organism record into database
        organism_id = organism_record.insert_into_database(cur)

        # connect organism record to compound record in database
        compound_record_id = compound_record.compound_id
        connect_compound_to_organism(cur, compound_record_id, organism_id)

        # insert biosynthetic class records into database
        for biosynthetic_class_record in biosynthetic_class_records:
            biosynthetic_class_id = biosynthetic_class_record.insert_into_database(cur)
            connect_compound_to_biosynthetic_class(cur, compound_record_id, biosynthetic_class_id)
