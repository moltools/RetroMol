#!/usr/bin/env python3

from retromol_database.duckdb import FINGERPRINT_SIZE
from retromol.model.rules import RuleSet
from retromol_fingerprint.fingerprint import TOKEN_UNK, Vocabulary, Fingerprinter
from retromol_alignment.scoring import create_tanimoto_scoring_matrix, HARDCODED_PK_SCORING
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import align
from retromol_alignment.pairwise import Converter


def main() -> None:
    monomers1 = [["acetic acid"],["C1", "PK"], ["C2", "PK"], ["B^S2^S", "PK"], ["B^S1", "PK"], ["D2^R", "PK"], ["D2^R", "PK"], ["B^R2^S", "PK"], ["C1", "PK"], ["B^R1", "PK"], ["B^S2^R", "PK"], ["C1", "PK"], ["C1", "PK"]]
    monomers2 = [["acetic acid"], ["C1", "PK"], ["C2", "PK"], ["B^S2^S", "PK"], ["B^S2^S", "PK"], ["D2", "PK"], ["C2", "PK"], ["B^R2^S", "PK"], ["C1", "PK"], ["B^R1", "PK"], ["B^S2^R", "PK"], ["B^S2^S", "PK"]]
    seq1 = [x[0] for x in monomers1]
    seq2 = [x[0] for x in monomers2]

    rules = RuleSet.load_default()
    vocab_tokens = set()
    for rule in rules.matching_rules:
        rule_name = rule.name
        vocab_tokens.add(rule_name)
        rule_pseudonyms = rule.pseudonyms
        vocab_tokens.update(set(rule_pseudonyms))
    vocab = Vocabulary(vocab_tokens)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)

    fp1 = fingerprinter.encode(monomers1)
    fp2 = fingerprinter.encode(monomers2)
    score = fingerprinter.similarity(fp1, fp2)
    print(score)

    records = []
    for rule in rules.matching_rules:
        rule_name = rule.name
        rule_smiles = rule.smiles
        records.append((rule_name, rule_smiles))
    scoring_matrix = create_tanimoto_scoring_matrix(
        records,
        radius=2,
        num_bits=2048,
        stereochemistry=True,
        self_score_tokens=[TOKEN_UNK, "A", "B", "C", "D"],
        self_score=1.0,
        hardcoded_scores=HARDCODED_PK_SCORING,
    )
    aligner = setup_aligner(scoring_matrix, mode="global")
    converter = Converter(
        to_identifier=lambda s: scoring_matrix.alphabet.index(s),
        from_identifier=lambda i: scoring_matrix.alphabet[i],
    )
    max_score1, _, _ = align(aligner, seq1, seq1, converter)
    max_score2, _, _ = align(aligner, seq2, seq2, converter)
    print(max_score1, max_score2)
    max_score = max(max_score1, max_score2)
    score, aln1, aln2 = align(aligner, seq1, seq2, converter)
    print(score/max_score)
    for a, b in zip(aln1, aln2):
        print(f"{a}\t{b}")



if __name__ == "__main__":
    main()
