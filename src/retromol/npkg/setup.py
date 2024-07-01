# -*- coding: utf-8 -*-

"""Setup NPKG constraints."""

import logging
import os

from tqdm import tqdm

from retromol.npkg.connection import Neo4jConnection
from retromol.npkg.helpers import validate_path
from retromol.npkg.nodes import (
    BioactivityLabel,
    BiosyntheticGeneCluster,
    Compound,
    Motif,
    MotifCode,
    Organism,
    Pathway,
    Protocluster,
)
from retromol.npkg.parsing import parse_compounds, parse_donphan, parse_mibig, parse_npatlas
from retromol.npkg.purge import purge_retrosynthesis_data


def set_constraints(conn: Neo4jConnection) -> None:
    """Set constraints on the Neo4j database.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :raises TypeError: If the connection is not a Neo4jConnection.
    """
    logger = logging.getLogger(__name__)
    logger.info("Setting constraints...")

    if not isinstance(conn, Neo4jConnection):
        msg = f"Expected a Neo4jConnection but received {type(conn)}"
        logger.error(msg)
        raise TypeError(msg)

    # Set node constraints.
    BioactivityLabel.set_constraints(conn)
    BiosyntheticGeneCluster.set_constraints(conn)
    Compound.set_constraints(conn)
    Motif.set_constraints(conn)
    MotifCode.set_constraints(conn)
    Organism.set_constraints(conn)
    Pathway.set_constraints(conn)
    Protocluster.set_constraints(conn)

    logger.info("Constraints set.")


def create_database(
    conn: Neo4jConnection,
    path_to_data: str,
    path_to_rxn: str,
    path_to_mon: str,
    set_constraints: bool = True,
    only_retrosynthesis: bool = False,
    num_workers: int = 1,
) -> None:
    """Create the NPKG database.

    :param conn: The Neo4j connection.
    :type conn: Neo4jConnection
    :param path_to_data: The path to the directory containing the data for the database construction.
    :type path_to_data: str
    :param path_to_rxn: The path to the JSON file containing the reactions.
    :type path_to_rxn: str
    :param path_to_mon: The path to the JSON file containing the monomer patterns.
    :type path_to_mon: str
    :param set_constraints: If True, set constraints on the database.
    :type set_constraints: bool
    :param only_retrosynthesis: Only reparse retrosynthesis data.
    :type only_retrosynthesis: bool
    :param num_workers: The number of workers to use for parsing.
    :type num_workers: int
    :raises TypeError: If the connection is not a Neo4jConnection.
    """
    logger = logging.getLogger(__name__)

    if not only_retrosynthesis:
        logger.info("Creating the NPKG database...")

        # Check if all necessary files are present in the data directory.
        path_donphan = os.path.join(path_to_data, "donphan.csv")
        validate_path(path_donphan)

        path_npatlas = os.path.join(path_to_data, "npatlas.json")
        validate_path(path_npatlas)

        path_mibig = os.path.join(path_to_data, "mibig")
        validate_path(path_mibig)

        # Check database connection.
        if not isinstance(conn, Neo4jConnection):
            msg = f"Expected a Neo4jConnection but received {type(conn)}"
            logger.error(msg)
            raise TypeError(msg)

        # Set constraints.
        if set_constraints:
            set_constraints(conn)

        # Load data.
        parse_npatlas(conn, path_npatlas)
        parse_donphan(conn, path_donphan)
        parse_mibig(conn, path_mibig)

        # Parse compounds.
        parse_compounds(conn, path_to_rxn, path_to_mon, num_workers)

        logger.info("NPKG database created.")

    else:
        logger.info("Re-parsing retrosynthesis data...")

        # Purge retrosynthesis data.
        purge_retrosynthesis_data(conn)

        # Set constraints.
        if set_constraints:
            set_constraints(conn)

        # Load data.
        parse_compounds(conn, path_to_rxn, path_to_mon, num_workers)

        logger.info("Retrosynthesis data re-parsed.")
