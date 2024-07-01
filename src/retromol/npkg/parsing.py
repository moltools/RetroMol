# -*- coding: utf-8 -*-

"""Parsing data for the NPKG database."""

import json
import logging
import multiprocessing as mp
import os
import typing as ty
from collections import defaultdict

from rdkit import Chem, RDLogger
from tqdm import tqdm

from retromol.npkg.connection import Neo4jConnection
from retromol.npkg.helpers import safe_get
from retromol.npkg.nodes import (
    BioactivityLabel,
    BiosyntheticGeneCluster,
    Compound,
    Motif,
    MotifCode,
    Organism,
    Pathway,
)
from retromol.retrosynthesis.chem import MolecularPattern, Molecule, ReactionRule
from retromol.retrosynthesis.helpers import timeout
from retromol.retrosynthesis.parsing import (
    parse_mol,
    parse_molecular_patterns,
    parse_reaction_rules,
)
from retromol.retrosynthesis.result import Result


@timeout(10)
def parse_mol_timed(
    record: ty.Tuple[str, str, Molecule, ty.List[ReactionRule], ty.List[MolecularPattern]]
) -> ty.Tuple[str, str, Result]:
    """Parse a molecule with a timeout.

    :param record: The record containing the compound_identifier, compound_source,
        molecule, reaction rules, and molecular patterns.
    :type record: Record
    :return: The result of the parsing.
    :rtype: Result
    """
    logger = logging.getLogger(__name__)

    RDLogger.DisableLog("rdApp.*")

    try:
        compound_identifier, compound_source, mol, reactions, monomers = record
        return compound_identifier, compound_source, parse_mol(mol, reactions, monomers)

    except Exception as e:
        msg = f"Failed to parse {mol.name} and raised {e.__class__.__name__}: {e}"
        logger.error(msg)
        return compound_identifier, compound_source, Result(mol.name, mol.compiled, False)


def parse_batch(conn: Neo4jConnection, num_workers: int, jobs: list) -> None:
    """Parse a batch of compounds with RetroMol.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param num_workers: The number of workers to use.
    :type num_workers: int
    :param jobs: The list of jobs to parse.
    :type jobs: list
    """
    logger = logging.getLogger(__name__)

    # Determine number of workers.
    num_workers = min(num_workers, mp.cpu_count())

    # Create pool of workers.
    with mp.Pool(num_workers) as pool:
        for compound_identifier, compound_source, result in tqdm(
            pool.imap_unordered(parse_mol_timed, jobs),
            total=len(jobs),
            desc="Parsing batch of compounds with RetroMol",
            leave=False,
        ):
            # Skip if result is not successful.
            if result.success is False:
                continue

            # Retrieve information from result.
            applied_reactions = result.applied_reactions
            sequences = result.sequences

            # Add biosynthetic fingerprint to database.
            for sequence in sequences:
                motif_code = sequence["motif_code"]

                if len(motif_code) < 3:
                    logger.warning(f"Skipping motif code with less than 3 motifs: {motif_code}")
                    continue

                # Create sub query for motif code.
                mc_query, mc_query_props = MotifCode.create_sub_query(
                    compound_identifier=compound_identifier,
                    compound_source=compound_source,
                    calculated=True,
                    applied_reactions=applied_reactions,
                    src=json.dumps(motif_code),
                )

                # Create full query.
                query_begin = f"CREATE {mc_query}-[:START]->"
                query_end = "-[:NEXT]->".join(
                    [
                        Motif.from_string(motif_index, motif)
                        for motif_index, motif in enumerate(motif_code)
                    ]
                )
                query_full = query_begin + query_end

                # Execute query.
                conn.query(query_full, mc_query_props)

                # Connect MotifCode to Compound node.
                conn.query(
                    (
                        "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "
                        "MATCH (mc:MotifCode {compound_identifier: $compound_identifier, compound_source: $compound_source}) "
                        "MERGE (c)-[:HAS_MOTIF_CODE]->(mc)"
                    ),
                    {
                        "compound_identifier": compound_identifier,
                        "compound_source": compound_source,
                    },
                )


def parse_compounds(
    conn: Neo4jConnection, path_to_rxn: str, path_to_mon: str, num_workers: int = 1
) -> None:
    """Parse compounds from the Neo4j database with RetroMol.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param path_to_rxn: The path to the reaction rules in json format.
    :type path_to_rxn: str
    :param path_to_mon: The path to the molecular patterns in json format.
    :type path_to_mon: str
    :param num_workers: The number of workers to use.
    :type num_workers: int
    """
    logger = logging.getLogger(__name__)

    # Disable RDKit logger.
    RDLogger.DisableLog("rdApp.*")

    # Parse reactions.
    logger.info("Parsing reaction rules...")
    reactions_src = json.load(open(path_to_rxn, "r", encoding="utf-8"))
    reactions = parse_reaction_rules(reactions_src)
    logger.info(f"Parsed {len(reactions)} reaction rules.")

    # Parse monomers.
    logger.info("Parsing molecular patterns...")
    monomers_src = json.load(open(path_to_mon, "r", encoding="utf-8"))
    monomers = parse_molecular_patterns(monomers_src)
    logger.info(f"Parsed {len(monomers)} molecular patterns.")

    # Retrieve compoounds from 'Polyketides' and 'Amino acids and Peptides' pathways.
    results = conn.query(
        (
            "MATCH (c:Compound)-[:BELONGS_TO]->(p:Pathway) "
            "WHERE p.name IN ['Polyketides', 'Amino acids and Peptides'] "
            "RETURN c.identifier AS identifier, c.source AS source, c.inchi AS inchi"
        ),
        batch_size=1,
    )

    # Collect jobs for multiprocessing.
    jobs = []

    # Unpack results and parse with RetroMol.
    for result in tqdm(results, desc="Parsing compounds with RetroMol", leave=False):

        compound_identifier = result["identifier"]
        compound_source = result["source"]
        compound_inchi = result["inchi"]

        # Parse compound with RDKit into RetroMol Molecule.
        try:
            mol = Chem.MolFromInchi(compound_inchi)

            if mol is None:
                raise Exception(f"Failed to parse InChI")

            smi = Chem.MolToSmiles(mol)
            rec = Molecule(compound_identifier, smi)

        except Exception as e:
            logger.error(f"Failed to parse compound {compound_identifier} with error: {e}")
            continue

        # Apply retrosynthesis to the compound record.
        job = (compound_identifier, compound_source, rec, reactions, monomers)
        jobs.append(job)

        # Collect jobs in batches of 1000.
        if len(jobs) == 1000:
            parse_batch(conn, num_workers, jobs)
            jobs = []

    # Parse remaining jobs.
    if jobs:
        parse_batch(conn, num_workers, jobs)


def parse_donphan(conn: Neo4jConnection, path: str) -> None:
    """Parse the DONPHAN data.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param path: The path to the DONPHAN data in csv format.
    :type path: str
    """
    # Disable RDKit logger.
    RDLogger.DisableLog("rdApp.*")

    # Collect bioactivity data.
    bioactivity_library = defaultdict(set)

    with open(path, "r", encoding="utf-8") as file_open:
        header = file_open.readline().strip().split(",")

        for line in tqdm(file_open, desc="Parsing DONPHAN data", leave=False):
            items = list(zip(header, line.strip().split(",")))

            smiles = items[1][1]
            inchikey = Chem.MolToInchiKey(Chem.MolFromSmiles(smiles))
            inchikey_connectivity = inchikey.split("-")[0]

            bioactivities = items[2:]
            bioactivities = set([h for h, b in bioactivities if b != ""])
            bioactivities = [x[:-3] if x.endswith("_np") else x for x in bioactivities]

            bioactivity_library[inchikey_connectivity].update(bioactivities)

    # Loop through the bioactivity library and add bioactivity nodes and relationships.
    for inchikey_connectivity, bioactivities in tqdm(
        bioactivity_library.items(), desc="Parsing DONPHAN data", leave=False
    ):  # noqa: E501

        # First retrieve compounds with the given connectivity.
        results = conn.query(
            (
                "MATCH (c:Compound {inchikey_connectivity: $inchikey_connectivity}) "
                "RETURN c.identifier AS identifier, c.source AS source"
            ),
            {"inchikey_connectivity": inchikey_connectivity},
        )

        # Annotate retrieved compounds with bioactivity data.
        for result in results:

            compound_identifier = result["identifier"]
            compound_source = result["source"]

            for bioactivity in bioactivities:
                bioactivity_source = "donphan"

                # Create bioactivity node.
                BioactivityLabel.create(conn, name=bioactivity, source=bioactivity_source)

                # Create relationship between compound and bioactivity.
                conn.query(
                    (
                        "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "  # noqa: E501
                        "MATCH (b:BioactivityLabel {name: $bioactivity, source: $bioactivity_source}) "  # noqa: E501
                        "MERGE (c)-[:HAS_BIOACTIVITY]->(b)"
                    ),
                    {
                        "compound_identifier": compound_identifier,
                        "compound_source": compound_source,
                        "bioactivity": bioactivity,
                        "bioactivity_source": bioactivity_source,
                    },
                )


def parse_mibig(conn: Neo4jConnection, path: str) -> None:
    """Parse the MIBIG data.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param path: The path to the MIBIG data in json format.
    :type path: str
    """
    paths = [
        os.path.join(path, file_path)
        for file_path in os.listdir(path)
        if file_path.endswith(".json")
    ]

    for path in tqdm(paths, desc="Parsing MIBIG data", leave=False):

        # Get cluster data.
        data = json.load(open(path, "r", encoding="utf-8"))
        cluster = data["cluster"]

        # Retrieve BGC information.
        bgc_identifier = cluster["mibig_accession"]
        bgc_source = "mibig"

        # Create BGC node.
        BiosyntheticGeneCluster.create(conn, identifier=bgc_identifier, source=bgc_source)

        # Parse produced compounds from the cluster, link BGC to compounds.
        if compounds := cluster.get("compounds", None):
            for compound in compounds:

                # Retrieve annotated bioactivities.
                bioactivities = []
                if chem_acts := compound.get("chem_acts", None):
                    for chem_act in chem_acts:
                        bioactivities.append(chem_act["activity"])

                # Retrieve compound NPAtlas identifier, if any, and use it to
                # link to the compound.
                if database_ids := compound.get("database_id", None):
                    for database_id in database_ids:

                        compound_source, compound_identifier = database_id.split(":")
                        if compound_source == "npatlas":

                            # Link compound to BGC.
                            conn.query(
                                (
                                    "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "  # noqa: E501
                                    "MATCH (b:BiosyntheticGeneCluster {identifier: $bgc_identifier, source: $bgc_source}) "  # noqa: E501
                                    "MERGE (b)-[:PRODUCES]->(c) "
                                    "MERGE (c)-[:PRODUCED_BY]->(b)"
                                ),
                                {
                                    "compound_identifier": compound_identifier,
                                    "compound_source": compound_source,
                                    "bgc_identifier": bgc_identifier,
                                    "bgc_source": bgc_source,
                                },
                            )

                            # Link compound to bioactivities.
                            for bioactivity in bioactivities:
                                bioactivity_source = "mibig"

                                # Create bioactivity node.
                                BioactivityLabel.create(
                                    conn, name=bioactivity, source=bioactivity_source
                                )  # noqa: E501

                                # Create relationship between compound and bioactivity.
                                conn.query(
                                    (
                                        "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "  # noqa: E501
                                        "MATCH (b:BioactivityLabel {name: $bioactivity, source: $bioactivity_source}) "  # noqa: E501
                                        "MERGE (c)-[:HAS_BIOACTIVITY]->(b)"
                                    ),
                                    {
                                        "compound_identifier": compound_identifier,
                                        "compound_source": compound_source,
                                        "bioactivity": bioactivity,
                                        "bioactivity_source": bioactivity_source,
                                    },
                                )


def parse_npatlas(conn: Neo4jConnection, path: str) -> None:
    """Parse the NPAtlas data.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param path: The path to the NPAtlas data in json format.
    :type path: str
    """
    with open(path, "r") as file:
        data = json.load(file)

    for record in tqdm(data, desc="Parsing NPAtlas data", leave=False):

        # Get producing organism information.
        organism = record["origin_organism"]
        organism_ncbi_id = organism["taxon"]["ncbi_id"]
        organism_type = organism["type"]
        organism_genus = organism["genus"]
        organism_species = organism["species"]

        # Get pathway information.
        pathways = safe_get(record, ["npclassifier", "pathway_results"], [])
        pathways_source = "npclassifier"

        # Get compound information.
        compound_source = "npatlas"
        compound_identifier = record["npaid"]
        compound_inchikey = record["inchikey"]
        compound_inchikey_connectivity = compound_inchikey.split("-")[0]
        compound_inchi = record["inchi"]

        # Create compound node.
        Compound.create(
            conn,
            identifier=compound_identifier,
            source=compound_source,
            inchikey=compound_inchikey,
            inchikey_connectivity=compound_inchikey_connectivity,
            inchi=compound_inchi,
        )

        # Create organism node, and relationship between compound and organism.
        Organism.create(
            conn,
            ncbi_id=organism_ncbi_id,
            type=organism_type,
            genus=organism_genus,
            species=organism_species,
        )

        conn.query(
            (
                "MATCH (c:Compound {identifier: $identifier, source: $source}) "
                "MATCH (o:Organism {ncbi_id: $ncbi_id}) "
                "MERGE (c)-[:PRODUCED_BY]->(o) "
                "MERGE (o)-[:PRODUCES]->(c)"
            ),
            {
                "identifier": compound_identifier,
                "source": compound_source,
                "ncbi_id": organism_ncbi_id,
            },
        )

        # Create pahtway node, and relationship between compound and pathway.
        for pathway in pathways:
            Pathway.create(
                conn,
                name=pathway,
                source=pathways_source,
            )

            conn.query(
                (
                    "MATCH (c:Compound {identifier: $compound_identifier, source: $compound_source}) "  # noqa: E501
                    "MATCH (p:Pathway {name: $pathway_name, source: $pathway_source}) "  # noqa: E501
                    "MERGE (c)-[:BELONGS_TO]->(p) "
                    "MERGE (p)-[:CONTAINS]->(c)"
                ),
                {
                    "compound_identifier": compound_identifier,
                    "compound_source": compound_source,
                    "pathway_name": pathway,
                    "pathway_source": pathways_source,
                },
            )
