# -*- coding: utf-8 -*-

"""This module contains functions for quering a query."""

import json
import os
import re
import typing as ty

import neo4j
import numpy as np
from flask import Blueprint, Response, request
from versalign.motif import Motif, Gap
from versalign.msa import multiple_sequence_alignment
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence
from tqdm import tqdm

from retromol.retrosynthesis.alignment import OtherMotif, PolyketideMotif, sequence_from_motif_string_list

from .common import NEO4J_PASSWORD, NEO4J_URI, NEO4J_USER, fail, success, warning


# SCORING_MATRIX = pd.read_csv(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data/scoring_matrix.csv"), index_col=0)
SCORING_MATRIX_LABELS = np.loadtxt(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data/scoring_matrix.csv"), delimiter=",", max_rows=1, dtype=str)
SCORING_MATRIX = np.loadtxt(os.path.join(os.path.dirname(os.path.dirname(__file__)), "data/scoring_matrix.csv"), delimiter=",", skiprows=1)


# TODO: This is a quick fix, bug in versalign.
def pad_results(motif_codes: ty.List[ty.Dict[str, ty.Any]]) -> ty.List[ty.Dict[str, ty.Any]]:
    """Make sure all aligned sequences are of equal length. If not, pad at end with gaps."""
    max_length = max([len(x["motifCode"]) for x in motif_codes])
    for item in motif_codes:
        if len(item["motifCode"]) < max_length:
            item["motifCode"] += [{"motifType": "gap", "polyketideType": "Any", "polyketideDecor": "Any", "peptideSource": "Any", "peptideCid": "Any"}] * (max_length - len(item["motifCode"]))
    return motif_codes


# TODO: This is a temporary solution.
def score_func(a: Motif, b: Motif) -> int:
    """Score function for pairwise alignment.

    :param a: Motif A.
    :type a: Motif
    :param b: Motif B.
    :type b: Motif
    :return: Score.
    :rtype: int

    TODO: Dynamically create score function based on supplied monomers.
    """
    if a == b:
        return 4
    
    elif isinstance(a, PolyketideMotif) and isinstance(b, PolyketideMotif):
        if (
            a.type == b.type
            and (
                a.decoration == b.decoration
                or a.decoration == "Any"
                or b.decoration == "Any"
            )
        ):
            return 3
        elif (
            a.type == b.type
            and (a.decoration == "Any" or b.decoration == "Any")
        ):
            return 3
        elif a.type == "Any" or b.type == "Any":
            return 2
        elif a.decoration == b.decoration:
            return 1
        else:
            0
    
    elif isinstance(a, OtherMotif) and isinstance(b, OtherMotif):

        if a.cid == "Any" and b.cid == "Any":
            return 1

        try:
            cid_a = str(a.cid)
            cid_b = str(b.cid)

            idx_cid_a = SCORING_MATRIX_LABELS.tolist().index(cid_a)
            idx_cid_b = SCORING_MATRIX_LABELS.tolist().index(cid_b)

            score = SCORING_MATRIX[idx_cid_a, idx_cid_b]
        
        except ValueError:
            return 0

        return score
    
    return -1


def compile_motif_query_item(motif: ty.Dict[str, ty.Any]) -> ty.List[ty.Dict[str, ty.Any]]:
    """Compile motif query item.

    :param motif: Motif to compile.
    :type motif: ty.Dict[str, ty.Any]
    :return: Compiled motif query item.
    :rtype: ty.List[ty.Dict[str, ty.Any]]
    """
    if motif["motifType"] == "polyketide":
        polyketide_type = motif["polyketideType"]
        polyketide_decor = motif["polyketideDecor"]

        # TODO: quick fix to remove stereochem for now
        polyketide_type = re.sub(r'[rR|sS]', '', str(polyketide_type))
        polyketide_decor = re.sub(r'[rR|sS]', '', str(polyketide_decor))
        polyketide_decor = int(polyketide_decor) if polyketide_decor != "Any" else "Any"

        return PolyketideMotif(
            type=polyketide_type,
            decoration=polyketide_decor
        )
    elif motif["motifType"] == "peptide":
        return OtherMotif(
            source=motif["peptideSource"],
            cid=motif["peptideCid"]
        )
    else:
        raise ValueError("Invalid motif type.")
    

def get_bioactivity_labels(session: neo4j.Session, compound_identifier: str) -> ty.List[str]:
    """Get bioactivity labels for a compound.
    
    :param session: Neo4j session.
    :type session: neo4j.Session
    :param compound_identifier: Compound identifier.
    :type compound_identifier: str
    :return: Bioactivity labels.
    :rtype: ty.List[str]
    """
    bioactivities = []
    query = "MATCH (c:Compound {identifier: $identifier})-[:HAS_BIOACTIVITY]->(b:BioactivityLabel) RETURN b"
    result = session.run(query, identifier=compound_identifier)
    for record in result:
        bioactivities.append(record["b"]["name"])
    return bioactivities


def get_genus_labels(session: neo4j.Session, compound_identifier: str) -> ty.List[str]:
    """Get genus labels for a compound.
    
    :param session: Neo4j session.
    :type session: neo4j.Session
    :param compound_identifier: Compound identifier.
    :type compound_identifier: str
    :return: Genus labels.
    :rtype: ty.List[str]
    """
    genus = []
    query = "MATCH (c:Compound {identifier: $identifier})-[:PRODUCED_BY]->(o:Organism) RETURN o"
    result = session.run(query, identifier=compound_identifier)
    for record in result:
        genus.append(record["o"]["genus"])
    return genus


def construct_statement_query_space(
    start_statement: str, 
    return_statement: str,
    query_against_molecules: bool,
    query_against_protoclusters: bool,
    selected_bioactivity_labels: ty.List[str],
    selected_organism_labels: ty.List[str]
) -> str:
    """Construct statement for defining the query space.

    :param start_statement: Start statement.
    :type start_statement: str
    :param return_statement: Return statement.
    :type return_statement: str
    :param query_against_molecules: Query against molecules.
    :type query_against_molecules: bool
    :param query_against_protoclusters: Query against protoclusters.
    :type query_against_protoclusters: bool
    :param selected_bioactivity_labels: Selected bioactivity labels.
    :type selected_bioactivity_labels: ty.List[str]
    :param selected_organism_labels: Selected organism labels.
    :type selected_organism_labels: ty.List[str]
    :return: Statement for defining the query space.
    :rtype: str
    """
    statement = start_statement

    if query_against_molecules and query_against_protoclusters:
        statement += "((b)<-[:HAS_MOTIF_CODE]-(:Compound) OR (b)<-[:HAS_MOTIF_CODE]-(:Protocluster)) "
    elif query_against_molecules:
        statement += "(b)<-[:HAS_MOTIF_CODE]-(:Compound) "
    elif query_against_protoclusters:
        statement += "(b)<-[:HAS_MOTIF_CODE]-(:Protocluster) "
    else:
        raise ValueError("No query space specified.") 
    
    # Only molecules are annotated with bioactivity labels.
    if selected_bioactivity_labels and query_against_molecules:
        for label in selected_bioactivity_labels:
            statement += f"AND (b)<-[:HAS_MOTIF_CODE]-(:Compound)-[:HAS_BIOACTIVITY]->(:BioactivityLabel {{name: '{label}'}}) "

    # Only molecules are annotated with organism labels.
    if selected_organism_labels and query_against_molecules:
        for ncbi_id in selected_organism_labels:
            statement += f"AND (b)<-[:HAS_MOTIF_CODE]-(:Compound)-[:PRODUCED_BY]->(:Organism {{ncbi_id: {ncbi_id}}}) "

    statement += return_statement

    return statement


def prep_query_results(
    session: neo4j.Session,
    results: ty.List[Sequence]
) -> ty.List[ty.Dict[str, ty.Any]]:
    """Prepare query results.

    :param session: Neo4j session.
    :type session: neo4j.Session
    :param results: Query results.
    :type results: ty.List[Sequence]
    :return: Prepared query results.
    :rtype: ty.List[ty.Dict[str, ty.Any]]
    """
    motif_codes = []

    for item in results:
        motif_code_identifier = item.identifier
        bioactivities = get_bioactivity_labels(session, motif_code_identifier)
        genus = get_genus_labels(session, motif_code_identifier)

        motif_code = []
        for motif in item._motifs:
            if isinstance(motif, PolyketideMotif):
                motif_code.append({
                    "motifType": "polyketide",
                    "polyketideType": motif.type,
                    "polyketideDecor": motif.decoration,
                    "peptideSource": "Any",
                    "peptideCid": "Any",
                })
            
            elif isinstance(motif, OtherMotif):
                motif_code.append({
                    "motifType": "peptide",
                    "polyketideType": "Any",
                    "polyketideDecor": "Any",
                    "peptideSource": motif.source,
                    "peptideCid": motif.cid,
                })

            elif isinstance(motif, Gap):
                motif_code.append({
                    "motifType": "gap",
                    "polyketideType": "Any",
                    "polyketideDecor": "Any",
                    "peptideSource": "Any",
                    "peptideCid": "Any",
                })

            else:
                print(motif.__class__, motif)
                raise ValueError("Invalid motif type retrieved from database.")
            
        motif_codes.append({
            "motifCode": motif_code,
            "motifCodeIdentifier": motif_code_identifier,
            "bioactivities": bioactivities,
            "genus": genus
        })

    return motif_codes


def pairwise_match(
    query: ty.List[ty.Dict[str, ty.Any]], 
    gap_penalty: int, 
    end_gap_penalty: int,
    max_num_matches: int,
    alignment_type: PairwiseAlignment,
    min_match_length: int,
    max_match_length: int,
    query_against_molecules: bool,
    query_against_protoclusters: bool,
    selected_bioactivity_labels: ty.List[str],
    selected_organism_labels: ty.List[str]
) -> ty.Dict[str, ty.Any]:
    """Pairwise match query against database.
    
    :param query: Query to match.
    :type query: ty.List[ty.Dict[str, ty.Any]]
    :param gap_penalty: Gap penalty.
    :type gap_penalty: int
    :param end_gap_penalty: End gap penalty.
    :type end_gap_penalty: int
    :param max_num_matches: Maximum number of matches.
    :type max_num_matches: int
    :param alignment_type: Alignment type.
    :type alignment_type: PairwiseAlignment
    :param min_match_length: Minimum match length.
    :type min_match_length: int
    :param max_match_length: Maximum match length.
    :type max_match_length: int
    :param query_against_molecules: Query against molecules.
    :type query_against_molecules: bool
    :param query_against_protoclusters: Query against protoclusters.
    :type query_against_protoclusters: bool
    :param selected_bioactivity_labels: Selected bioactivity labels.
    :type selected_bioactivity_labels: ty.List[str]
    :param selected_organism_labels: Selected organism labels.
    :type selected_organism_labels: ty.List[str]
    :return: Match results.
    :rtype: ty.Dict[str, ty.Any]
    """
    # Connect to database.
    if NEO4J_USER and NEO4J_PASSWORD:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    else:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI)

    # Recompile query for matching. Check if query is non-ambiguous.
    compiled_query = []
    for motif in query:
        if len(motif) != 1:
            raise ValueError("Query is ambiguous. Motif must have a single identity.")
        motif = motif[0]
        compiled_query.append(compile_motif_query_item(motif))

    compiled_query = Sequence("Query", compiled_query)

    # Match query against database.
    top_scores, top_seqs = [], []
    with driver.session() as session:

        # Constructy statement for defining the query space.
        statement = construct_statement_query_space(
            "MATCH (b:MotifCode) WHERE ",
            "RETURN b",
            query_against_molecules=query_against_molecules,
            query_against_protoclusters=query_against_protoclusters,
            selected_bioactivity_labels=selected_bioactivity_labels,
            selected_organism_labels=selected_organism_labels
        )
        result = session.run(statement, fetch_size=1)

        for record in tqdm(result):
            name = record["b"]["compound_identifier"]
            motif_code = record["b"]["src"]
            motif_code = json.loads(motif_code)

            # Filter on min/max match length.
            if (
                len(motif_code) < min_match_length 
                or len(motif_code) > max_match_length
            ):
                continue
            
            # Parse motif code.
            try:
                sequence = sequence_from_motif_string_list(name, motif_code)
            except Exception as e:
                print(e) # Put this in because I had some erroneous motif codes in my local db
                continue
            
            # Define alignment options.
            if alignment_type == PairwiseAlignment.SMITH_WATERMAN:
                options = {"gap_penalty": gap_penalty}
            elif alignment_type == PairwiseAlignment.NEEDLEMAN_WUNSCH:
                options = {"gap_penalty": gap_penalty, "end_gap_penalty": end_gap_penalty}
            else:
                raise ValueError("Unsupported alignment strategy.")

            # Perform pairwise alignment.
            _, _, score = align_pairwise(
                compiled_query, 
                sequence, 
                score_func, 
                alignment_type, 
                options
            )

            # Update top scores.
            if len(top_scores) < max_num_matches:
                top_scores.append(score)
                top_seqs.append(sequence)

            elif score > min(top_scores):
                index = top_scores.index(min(top_scores))
                top_scores[index] = score
                top_seqs[index] = sequence

            else:
                continue

    # Sort and discard scores.
    top_items = sorted(zip(top_seqs, top_scores), key=lambda x: x[1], reverse=True)
    top_items = [x[0] for x in top_items]  # Discard scores.

    # Perform MSA.
    top_items.append(compiled_query)
    msa = multiple_sequence_alignment(top_items, gap_penalty, end_gap_penalty, score_func)

    # Prepare query results.
    with driver.session() as session:
        motif_codes = prep_query_results(session, msa)

    motif_codes = pad_results(motif_codes)
    return {"motif_codes": motif_codes}


def pattern_match(
    query: ty.List[ty.Dict[str, ty.Any]], 
    gap_penalty: int,
    end_gap_penalty: int,
    max_num_matches: int,
    min_match_length: int,
    max_match_length: int,
    query_against_molecules: bool,
    query_against_protoclusters: bool,
    selected_bioactivity_labels: ty.List[str],
    selected_organism_labels: ty.List[str],
    query_has_leading_modules: bool,
    query_has_trailing_modules: bool
) -> ty.Dict[str, ty.Any]:
    """Pairwise match query against database.
    
    :param query: Query to match.
    :type query: ty.List[ty.Dict[str, ty.Any]]
    :param gap_penalty: Gap penalty.
    :type gap_penalty: int
    :param end_gap_penalty: End gap penalty.
    :type end_gap_penalty: int
    :param max_num_matches: Maximum number of matches.
    :type max_num_matches: int
    :param min_match_length: Minimum match length.
    :type min_match_length: int
    :param max_match_length: Maximum match length.
    :type max_match_length: int
    :param query_against_molecules: Query against molecules.
    :type query_against_molecules: bool
    :param query_against_protoclusters: Query against protoclusters.
    :type query_against_protoclusters: bool
    :param selected_bioactivity_labels: Selected bioactivity labels.
    :type selected_bioactivity_labels: ty.List[str]
    :param selected_organism_labels: Selected organism labels.
    :type selected_organism_labels: ty.List[str]
    :param query_has_leading_modules: Query has leading modules.
    :type query_has_leading_modules: bool
    :param query_has_trailing_modules: Query has trailing modules.
    :type query_has_trailing_modules: bool
    :return: Match results.
    :rtype: ty.Dict[str, ty.Any]
    """
    # Connect to database.
    if NEO4J_USER and NEO4J_PASSWORD:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))
    else:
        driver = neo4j.GraphDatabase.driver(NEO4J_URI)
    
    def compile_path_query(num_match_items: int, has_leading_modules: bool) -> str:
        """Compile the path query.
        
        :param num_match_items: The number of match items.
        :type num_match_items: int
        :param has_leading_modules: Whether the sequence has leading modules.
        :type has_leading_modules: bool
        :return: The compiled path query.
        :rtype: str
        """
        query = []

        if has_leading_modules:
            module_range = range(2, num_match_items + 1)
        else:
            module_range = range(1, num_match_items)

        for i in module_range:
            query.append(f"(u{i})-[:NEXT]->")

        if has_leading_modules:
            query.append(f"(u{num_match_items + 1})")
        else:
            query.append(f"(u{num_match_items})")

        query = "".join(query)

        return query

    def query_item_to_addendum(index: int, query_item: ty.Dict[str, ty.Any]) -> str:
        """Convert a query item to a statement addendum.
        
        :param index: The index of the match item.
        :type index: int
        :param query_item: The match item.
        :type query_item: ty.Dict[str, ty.Any]
        :return: The query addendum.
        :rtype: str
        """
        substatement = []
        for option in query_item:

            if option["motifType"] == "polyketide":
                polyketide_type = option["polyketideType"] # A-D or Any
                polyketide_decor = option["polyketideDecor"] # 1-12 or Any

                # TODO: quick fix to remove stereochem for now
                polyketide_type = re.sub(r'[rR|sS]', '', str(polyketide_type))
                polyketide_decor = re.sub(r'[rR|sS]', '', str(polyketide_decor))
                polyketide_decor = int(polyketide_decor) if polyketide_decor != "Any" else "Any"

                if polyketide_type == "Any":
                    polyketide_type = None
                
                if polyketide_decor == "Any":
                    polyketide_decor = None

                subquery = \
                    (f"(u{index}.motif_type = 'polyketide'") + \
                    (f" AND u{index}.polyketide_type = '{polyketide_type}'" if polyketide_type is not None else "") + \
                    (f" AND u{index}.polyketide_decoration_type = {polyketide_decor})" if polyketide_decor is not None else ")")
                
                substatement.append(subquery)

            elif option["motifType"] == "peptide":
                peptide_source = option["peptideSource"]  # Ignored for now.
                peptide_cid = option["peptideCid"] # CID or Any

                if peptide_cid == "Any":
                    peptide_cid = None

                subquery = \
                    (f" (u{index}.motif_type = 'peptide'") + \
                    (f" AND u{index}.peptide_cid = '{peptide_cid}')" if peptide_cid is not None else ")")
                
                substatement.append(subquery)

            else:
                pass

        subquery = " OR ".join(substatement)
        
        return " AND (" + subquery + ")" if subquery != "" else ""

    if query_has_leading_modules:
        statement = (
            "MATCH (b:MotifCode)-[:START]->(u1)"
            " MATCH seq = (u1)-[:NEXT*]->(seqEnd)"
            " MATCH path = (u1)-[:NEXT*]->" + compile_path_query(len(query), True)
        )
    else:
        statement = (
            "MATCH (b:MotifCode)-[:START]->(u1)"
            " MATCH seq = (u1)-[:NEXT*]->(seqEnd)"
            " MATCH path = " + compile_path_query(len(query), False)
        )

    statement = construct_statement_query_space(
        statement + " WHERE ",
        "",
        query_against_molecules=query_against_molecules,
        query_against_protoclusters=query_against_protoclusters,
        selected_bioactivity_labels=selected_bioactivity_labels,
        selected_organism_labels=selected_organism_labels
    )

    if not query_has_leading_modules and not query_has_trailing_modules:
        statement += f" AND (NOT ()-[:NEXT]->(u1) AND NOT (u{len(query)})-[:NEXT]->())"
    elif not query_has_leading_modules:
        statement += " AND NOT ()-[:NEXT]->(u1)"
    elif not query_has_trailing_modules:
        statement += f" AND NOT (u{len(query) + 1})-[:NEXT]->()"

    for i, match_item in enumerate(query):
        if not query_has_leading_modules:
            index = i + 1
        else:
            index = i + 2

        addendum = query_item_to_addendum(index, match_item)
        if addendum != "":
            statement += addendum

    statement += (
        f" AND (length(seq) >= {min_match_length - 1}"
        f" AND length(seq) <= {max_match_length - 1}"
        f" AND NOT (seqEnd)-[:NEXT]->())"
    )

    statement += f" RETURN DISTINCT b LIMIT {max_num_matches}"

    with driver.session() as session:
        result = session.run(statement, fetch_size=1)

        sequences = []
        for record in result:
            name = record["b"]["compound_identifier"]
            motif_code = record["b"]["src"]
            motif_code = json.loads(motif_code)
            try:
                sequence = sequence_from_motif_string_list(name, motif_code)
            except Exception as e:
                print(e) # Put this in because I had some erroneous motif codes in my local db
                continue
            sequences.append(sequence)

        if len(sequences) == 0:
            return {"motif_codes": []}

        msa = multiple_sequence_alignment(sequences, gap_penalty, end_gap_penalty, score_func)
        motif_codes = prep_query_results(session, msa)
    
    motif_codes = pad_results(motif_codes)
    return {"motif_codes": motif_codes}


blueprint_query_submission = Blueprint("query_submission", __name__)
@blueprint_query_submission.route("/api/query_submission", methods=["POST"])
def query_submission() -> Response:
    """Query submission data."""
    
    # Unpack data.
    data = request.json
    
    try:
        query = data["query"]
        query_type = data["queryType"]
        alignment_type = data["alignmentType"]
        end_gap_penalty = data["endGapPenalty"]
        gap_penalty = data["gapPenalty"]
        min_match_length = data["minMatchLength"]
        max_match_length = data["maxMatchLength"]
        max_num_matches = data["maxNumMatches"]
        query_against_molecules = data["queryAgainstMolecules"]
        query_against_protoclusters = data["queryAgainstProtoclusters"]
        query_has_leading_modules = data["queryHasLeadingModules"]
        query_has_trailing_modules = data["queryHasTrailingModules"]
        selected_bioactivity_labels = data["selectedBioactivityLabels"]
        selected_organism_labels = data["selectedOrganismLabels"]

        assert len(query) > 0, "Input query is empty."
        assert query_type in ["match", "query"], "Invalid query type."
        assert min_match_length > 0, "Invalid minimum match length."
        assert max_match_length > 0, "Invalid maximum match length."
        assert max_match_length >= min_match_length, "Invalid match length range."
        assert max_num_matches > 0, "Invalid maximum number of matches."
        assert max_num_matches <= 100, "Maximum number of matches is too high."

        if alignment_type == "local":
            alignment_type = PairwiseAlignment.SMITH_WATERMAN
        elif alignment_type == "global":
            alignment_type = PairwiseAlignment.NEEDLEMAN_WUNSCH
        else:
            raise ValueError("Invalid alignment strategy.")
        
        selected_bioactivity_labels = [x["value"] for x in selected_bioactivity_labels]
        selected_organism_labels = [x["value"] for x in selected_organism_labels] # list of lists
        selected_organism_labels = [item for sublist in selected_organism_labels for item in sublist] # unpack list of lists

    except Exception as e:
        error_type = e.__class__.__name__
        return fail(f"Invalid input data ({error_type}): {e}")
    
    try:
        if query_type == "match":
            results = pairwise_match(
                query=query,
                gap_penalty=gap_penalty,
                end_gap_penalty=end_gap_penalty,
                max_num_matches=max_num_matches,
                alignment_type=alignment_type,
                min_match_length=min_match_length,
                max_match_length=max_match_length,
                query_against_molecules=query_against_molecules,
                query_against_protoclusters=query_against_protoclusters,
                selected_bioactivity_labels=selected_bioactivity_labels,
                selected_organism_labels=selected_organism_labels,
            )
            return success("Query ran successfully.", results)
        
        elif query_type == "query":
            results = pattern_match(
                query=query,
                gap_penalty=gap_penalty,
                end_gap_penalty=end_gap_penalty,
                max_num_matches=max_num_matches,
                min_match_length=min_match_length,
                max_match_length=max_match_length,
                query_against_molecules=query_against_molecules,
                query_against_protoclusters=query_against_protoclusters,
                selected_bioactivity_labels=selected_bioactivity_labels,
                selected_organism_labels=selected_organism_labels,
                query_has_leading_modules=query_has_leading_modules,
                query_has_trailing_modules=query_has_trailing_modules
            )
            return success("Query ran successfully.", results)
        
    except Exception as e:
        error_type = e.__class__.__name__
        return fail(f"Query failed ({error_type}): {e}")
    
    return fail("Endpoint not implemented yet!")
