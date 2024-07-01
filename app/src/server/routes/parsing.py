# -*- coding: utf-8 -*-

"""This module contains functions for parsing input data."""

import json
import os
import re
import typing as ty

import neo4j
from flask import Blueprint, Response, request
from rdkit import Chem

from retromol.npkg.antismash import get_antismash_data, parse_antismash_json
from retromol.retrosynthesis.parsing import (
    Molecule,
    parse_mol, 
    parse_molecular_patterns, 
    parse_reaction_rules
)

from .common import fail, warning, success


try:
    absolute_path = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    monomers = json.load(open(os.path.join(absolute_path, "data/monomers.json"), "r", encoding="utf-8"))
    reactions = json.load(open(os.path.join(absolute_path, "data/reactions.json"), "r", encoding="utf-8"))
    REACTIONS = parse_reaction_rules(reactions)
    MONOMERS = parse_molecular_patterns(monomers)
except Exception:
    REACTIONS = []
    MONOMERS = []


def motif_code_to_query(motif_code: ty.List[str]) -> ty.List[ty.List[str]]:
    """Convert motif code to a Cypher query.
    
    :param motif_code: List of motifs.
    :type motif_code: list[str]
    :return: List of queries.
    :rtype: list[list[str]]
    """
    wrapped_motif_code = []

    for motif in motif_code:
        if match := re.match(r"polyketide\|([A-D]|\*)(\d{1,2}|\*)", motif):

            if match.group(1) == "*":
                polyketide_type = "Any"
            else:
                polyketide_type = match.group(1)

            if match.group(2) == "*":
                polyketide_decor = "Any"
            else:
                polyketide_decor = int(match.group(2))
            
            parsed_motif = {
                "motifType": "polyketide",
                "polyketideType": polyketide_type,
                "polyketideDecor": polyketide_decor,
                "peptideSource": "Any",
                "peptideCid": "Any"
            }

        elif match := re.match(r"peptide\|(\w+)\|(.+)", motif):
            peptide_source = match.group(1)

            if match.group(2) == "*":
                peptide_cid = "Any"
            else:
                peptide_cid = match.group(2)

            parsed_motif = {
                "motifType": "peptide",
                "polyketideType": "Any",
                "polyketideDecor": "Any",
                "peptideSource": peptide_source,
                "peptideCid": peptide_cid
            }

        else:
            raise ValueError(f"Invalid motif: {motif}")    

        wrapped_motif_code.append([parsed_motif])  # Wrap individual motifs in a list.
    
    return wrapped_motif_code


blueprint_parse_submission = Blueprint("parse_submission", __name__)
@blueprint_parse_submission.route("/api/parse_submission", methods=["POST"])
def parse_submission() -> Response:
    """Parse submission data."""
    
    # Unpack data.
    data = request.json

    try:
        input_type = data["inputType"]
        input_value = data["inputValue"]

        assert input_type in ["smiles", "jobId", "json"], "Invalid input type."
        assert len(input_value) > 0, "Input value is empty."

    except Exception as e:
        error_type = e.__class__.__name__
        return fail(f"Invalid input data ({error_type}): {e}")
    
    # Parse data.
    try:

        ########################################################################
        #
        # Parse SMILES
        #
        ########################################################################
        if input_type == "smiles":
            molecule = Molecule("input", input_value)
            result = parse_mol(molecule, REACTIONS, MONOMERS)
            result.applied_reactions.sort()
            queries = [
                {   
                    "title": f"Retrosynthesis",
                    "queryType": "retrosynthesis",
                    "query": motif_code_to_query(seq["motif_code"]),
                    "metaData": {
                        "inputSmiles": Chem.MolToSmiles(result.mol),
                        "appliedRules": result.applied_reactions,
                        "mapping": seq["meta_data"]
                    }
                }
                for seq_index, seq in enumerate(result.sequences)
            ]

            if result.success is True:
                return success("Molecule parsed successfully!", {"queries": queries})
            else:
                return fail("Failed to parse molecule!")
        
        ########################################################################
        #
        # Parse job ID
        #
        ########################################################################
        elif input_type == "jobId":
            data = get_antismash_data(input_value)
            queries = parse_antismash_json(data)
            return success("Job parsed successfully!", {"queries": queries})
        
        ########################################################################
        #
        # Parse JSON
        #
        ########################################################################
        elif input_type == "json":
            data = json.loads(input_value)
            queries = parse_antismash_json(data)
            return success("JSON parsed successfully!", {"queries": queries})
    
    except Exception as e:
        error_type = e.__class__.__name__
        return fail(f"Failed to parse input ({error_type}): {e}")

    return warning("Endpoint not implemented yet!")
