import argparse
import json
import numpy as np
import os
from tqdm import tqdm

from retromol.npkg.connection import Neo4jConnection
from versalign.motif import Motif
from versalign.pairwise import PairwiseAlignment, align_pairwise
from versalign.sequence import Sequence


class A(Motif):
    def __eq__(self, other): return isinstance(other, A)
    def __str__(self): return "A"

class B(Motif):
    def __eq__(self, other): return isinstance(other, B)
    def __str__(self): return "B"

class C(Motif):
    def __eq__(self, other): return isinstance(other, C)
    def __str__(self): return "C"

class D(Motif):
    def __eq__(self, other): return isinstance(other, D)
    def __str__(self): return "D"

class Pep(Motif):
    def __eq__(self, other): return isinstance(other, Pep)
    def __str__(self): return "Pep"


def score_func(a, b):
    if a == b:
        return 3
    return -1


def cli():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default=".", help="The output directory.")
    parser.add_argument("--uri", type=str, default="bolt://localhost:7687", help="The URI for the Neo4j database.")
    parser.add_argument("--usr", type=str, default="neo4j", help="The user for the Neo4j database.")
    parser.add_argument("--pwd", type=str, default="password", help="The password for the Neo4j database.")
    return parser.parse_args()


def parse_sequence(identifier, seq_src):
    data = json.loads(seq_src)
    seq = []
    for x in data:
        if x.startswith("peptide"): seq.append(Pep())
        elif x.startswith("polyketide"):
            _, motif_code = x.split("|")
            if motif_code.startswith("A"): seq.append(A())
            elif motif_code.startswith("B"): seq.append(B())
            elif motif_code.startswith("C"): seq.append(C())
            elif motif_code.startswith("D"): seq.append(D())
            else: raise ValueError(f"Unknown motif code: {motif_code}")
        else: raise ValueError(f"Unknown motif code: {x}")
    return Sequence(identifier, seq)


def main():
    args = cli()
    conn = Neo4jConnection(args.uri, args.usr, args.pwd)
    labels, identifiers, seqs = [], [], []
    result = conn.query("MATCH (n:MotifCode) RETURN n")
    for record in result:
        label = record["n"]["compound_source"]  # npatlas / asdb4
        identifier = record["n"]["compound_identifier"]
        seq_src = record["n"]["src"]  # jsonified str sequence
        labels.append(label)
        identifiers.append(identifier)
        seqs.append(parse_sequence(identifier, seq_src))

    # Separate labels, identifiers, seqs based on label.
    npatlas_labels, npatlas_identifiers, npatlas_seqs = [], [], []
    asdb4_labels, asdb4_identifiers, asdb4_seqs = [], [], []

    for label, identifier, seq in zip(labels, identifiers, seqs):
        if label == "npatlas":
            npatlas_labels.append(label)
            npatlas_identifiers.append(identifier)
            npatlas_seqs.append(seq)
        elif label == "asdb4":
            asdb4_labels.append(label)
            asdb4_identifiers.append(identifier)
            asdb4_seqs.append(seq)
        else:
            raise ValueError(f"Unknown label: {label}")
        

    # target = Sequence("target", [B(), B(), A(), D(), B(), B()])

    out_file_path = os.path.join(args.out, "deorphanized.txt")
    with open(out_file_path, "w") as out_file:

        for target in tqdm(asdb4_seqs):
            num_peptides = sum([isinstance(x, Pep) for x in target._motifs])
            if num_peptides > 2:
                continue

            if len(target) < 6:
                continue

            max_pairwise_score = 3.0
            max_score = len(target) * max_pairwise_score

            top_scores = []
            for i, query in enumerate(npatlas_seqs):
                _, _, score = align_pairwise(target, query, score_func, PairwiseAlignment.SMITH_WATERMAN)
                score /= max_score
                # print(len(target), len(query), score)
                if (
                    score > ((len(target) - 2) * max_pairwise_score) / max_score
                    and len(query) > (len(target) - 2)
                    and len(query) < (len(target) + 2)
                ):
                    top_scores.append((query, score))

            top_scores = sorted(top_scores, key=lambda x: x[1], reverse=True)

            if len(top_scores) == 0:
                continue

            out_file.write(f"TARGET {target.identifier}: {target}\n")
            out_file.write(f"NUM_HITS: {len(top_scores)}\n")
            for k, (top_score, score) in enumerate(top_scores):
                k += 1
                out_file.write(f"HIT_{k}: {top_score.identifier} (score={round(score, 5)}): {top_score}\n")
            out_file.write("\n")

            # make sure written stuff is updated in file
            out_file.flush()


if __name__ == "__main__":
    main()
