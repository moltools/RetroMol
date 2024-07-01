# -*- coding: utf-8 -*-

"""This module contains functions for parsing input data."""

import json
import os
import re
import typing as ty
from collections import defaultdict

import neo4j
from flask import Blueprint, Response, request

from .common import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, fail, warning, success


blueprint_fetch_bioactivity_labels = Blueprint("fetch_bioactivity_labels", __name__)
@blueprint_fetch_bioactivity_labels.route("/api/fetch_bioactivity_labels", methods=["GET"])
def fetch_bioactivity_labels() -> Response:
    """Fetch bioactivity labels."""
    if NEO4J_USER and NEO4J_PASSWORD:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    else:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI)

    # Retrieve all bioactivity labels.
    labels = []
    with driver.session() as session:
        query = """
        MATCH (b:BioactivityLabel)
        RETURN b
        """
        result = session.run(query)
        for record in result:
            labels.append(record["b"]["name"])

    labels = list(set(labels))
    labels.sort()
    labels = [{"label": label, "value": label} for label in labels]
    payload = {"bioactivityLabels": labels}

    return success("Bioactivity labels fetched successfully!", payload)


blueprint_fetch_organism_labels = Blueprint("fetch_organism_labels", __name__)
@blueprint_fetch_organism_labels.route("/api/fetch_organism_labels", methods=["GET"])
def fetch_organism_labels() -> Response:
    """Fetch organism labels."""
    if NEO4J_USER and NEO4J_PASSWORD:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    else:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI)

    # Retrieve all producing organisms
    # retrieved_ncbi_ids = set()
    # producing_organisms = []
    producing_genus = defaultdict(set)
    with driver.session() as session:
        query = """
        MATCH (o:Organism)
        RETURN o
        """
        result = session.run(query)
        for record in result:
            genus = record["o"]["genus"]
            # species = record["o"]["species"]
            ncbi = record["o"]["ncbi_id"]
            # name = f"{genus} {species}"
            # if ncbi in retrieved_ncbi_ids:
            #     continue
            # producing_organisms.append({"label": name, "value": ncbi})
            # retrieved_ncbi_ids.add(ncbi)
            producing_genus[genus].add(ncbi)

    # producing_organisms.sort(key=lambda x: x["label"])

    producing_genus = [
        {"label": genus, "value": list(ncbi_ids)}
        for genus, ncbi_ids in producing_genus.items()
    ]
    producing_genus.sort(key=lambda x: x["label"])

    # payload = {"organismLabels": producing_organisms}
    payload = {"organismLabels": producing_genus}

    return success("Organism labels fetched successfully!", payload)
