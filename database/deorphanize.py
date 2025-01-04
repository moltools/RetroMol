#!/usr/bin/env python3
import argparse 
import os

from tqdm import tqdm
from versalign.motif import Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence

from populate_database import connect_to_database
from retrieval import SequenceMotif, parse_scoring_matrix

"""
TODO:
RUN DEORPHANIZATION AND GET ALIGNMENT SEQUENCES OUT AND NUM OF EXACT MATCHES!!!!
"""

def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=str, required=True, help="outdir deorphanization results")
    parser.add_argument("--scoring-matrix", type=str, required=True, help="scoring matrix for alignment")
    return parser.parse_args()


def get_compound_seqs(cur):
    query = """
    SELECT 
        c.compound_id, 
        STRING_AGG(m.motif_name, '|' ORDER BY psm.position) AS primary_sequence, 
        o.type AS organism_type, 
        o.genus AS organism_genus 
    FROM compounds c 
    JOIN compounds_primary_sequences cps ON c.compound_id = cps.compound_id 
    JOIN primary_sequences_motifs psm ON cps.primary_sequence_id = psm.primary_sequence_id 
    JOIN motifs m ON psm.motif_id = m.motif_name JOIN compounds_organisms co ON c.compound_id = co.compound_id 
    JOIN organisms o ON co.organism_id = o.id 
    GROUP BY c.compound_id, cps.primary_sequence_id, o.type, o.genus;
    """
    cur.execute(query)
    return cur.fetchall()


def get_protocluster_seqs(cur):
    query = """
    SELECT 
        p.id AS protocluster_id, 
        p.input_file AS input_file,
        STRING_AGG(m.motif_name, '|' ORDER BY psm.position) AS primary_sequence, 
        o.type AS organism_type, 
        o.genus AS organism_genus 
    FROM protoclusters p 
    JOIN protoclusters_primary_sequences pps ON p.id = pps.protocluster_id 
    JOIN primary_sequences_motifs psm ON pps.primary_sequence_id = psm.primary_sequence_id 
    JOIN motifs m ON psm.motif_id = m.motif_name 
    JOIN protoclusters_organisms po ON p.id = po.protocluster_id 
    JOIN organisms o ON po.organism_id = o.id
    GROUP BY p.id, pps.primary_sequence_id, o.type, o.genus;
    """
    cur.execute(query)
    return cur.fetchall()


def seq_repr(seq):
    return "|".join([str(m) for m in seq._motifs])


def main():
    args = cli()
    cur, conn = connect_to_database()
    scoring_matrix = parse_scoring_matrix(args.scoring_matrix)
    def score_func(a, b) -> int: return scoring_matrix[(str(a.name), str(b.name))]

    # get all protoclusters, and parse into alignable objects
    protocluster_seqs = get_protocluster_seqs(cur)
    protoclusters = []
    for protocluster_id, protocluster_input_file, primary_sequence, organism_type, organism_genus in tqdm(protocluster_seqs):
        primary_sequence = primary_sequence.split("|")
        protocluster = Sequence(protocluster_id, [SequenceMotif(m) for m in primary_sequence])
        protoclusters.append((protocluster_id, protocluster_input_file, protocluster, organism_type, organism_genus))
    print(f"Retrieved {len(protoclusters)} protoclusters")

    # retrieve compounds to deorphanize
    compound_seqs = get_compound_seqs(cur)
    for compound_id, primary_sequence, organism_type, organism_genus in tqdm(compound_seqs):
        primary_sequence = primary_sequence.split("|")
        
        # filter on primary sequence length
        if len(primary_sequence) < 7:
            continue

        # filter out protoclusters with same organism type
        targets = []
        for protocluster_id, protocluster_input_file, protocluster, protocluster_type, protocluster_genus in protoclusters:
            if (
                (protocluster_type == organism_type and protocluster_genus == organism_genus)
                or (organism_type == "unknown" and organism_genus == "unknown")
                or (protocluster_type == organism_type and organism_genus == "unknown")
            ):
                targets.append((protocluster_id, protocluster_input_file, protocluster, protocluster_type, protocluster_genus))
        if len(targets) == 0:
            continue

        # compile item for querying
        query = Sequence(compound_id, [SequenceMotif(m) for m in primary_sequence])

        # calcualte max score for query with score func
        max_score = sum([score_func(m, m) for m in query._motifs])

        # TODO: print correct protocluster id in out file
        # TODO: make sure sequences are within "" so they are not split by csv reader
        # TODO: add organism type and genus to out file
        # TODO: add raw score and max score to out file

        # align
        scores = []
        for target_id, target_input_file, target, target_type, target_genus in tqdm(targets, leave=False):
            aligned_a, aligned_b, score = align_pairwise(
                seq_a=query, 
                seq_b=target,
                score_func=score_func,
                algorithm=PairwiseAlignment.SMITH_WATERMAN,
                options={"gap_penalty": 3}
            )
            score /= max_score
            score = round(score, 4)
            scores.append((score, target_id, target_input_file, target, target_type, target_genus))
        scores = sorted(scores, key=lambda x: x[0], reverse=True)

        # create out file path with compound_id in it
        out_file = os.path.join(args.outdir, f"{compound_id}.csv")
        
        # write out top 50 results
        with open(out_file, "w") as f:
            f.write("compound_id,protocluster_id,protocluster_input_file,score,query_organism_type,query_organism_genus,target_organism_type,target_organism_genus,alignment_query,alignment_target\n")
            for score, target_id, target_input_file, target, target_type, target_genus in scores[:50]:
                f.write(f"{query.identifier},{target.identifier},{target_input_file},{score},{organism_type},{organism_genus},{target_type},{target_genus},\"{seq_repr(query)}\",\"{seq_repr(target)}\"\n")

    cur.close()
    conn.close()


if __name__ == "__main__":
    main()
