import argparse
from pathlib import Path

from tqdm import tqdm

from retromol_database.duckdb import FINGERPRINT_SIZE, open_database
from retromol.model.rules import RuleSet
from retromol_fingerprint.fingerprint import TOKEN_UNK, Vocabulary, Fingerprinter
from retromol_alignment.scoring import create_tanimoto_scoring_matrix
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import Converter
from retromol_alignment.ranking import rerank


PK_A_SUBTYPES = [
    "A1",
    "A2", "A2^R", "A2^S",
    "A3",
    "A4", "A4^R", "A4^S",
    "A5", "A5^R", "A5^S",
    "A6", "A6^R", "A6^S",
    "A7", "A7^R", "A7^S",
    "A8",
    "A9", "A9^R", "A9^S",
    "A10", "A10^R", "A10^S",
    "A11",
]
PK_B_SUBTYPES = [
    "B1", "B^R1", "B^S1",
    "B2", "B^R2^R", "B^R2^S", "B^S2^R", "B^S2^S", "B^R2", "B^S2",
    "B3", "B^R3", "B^S3",
    "B4", "B^R4^R", "B^R4^S", "B^S4^R", "B^S4^S",
    "B5", "B^R5^R", "B^R5^S", "B^S5^R", "B^S5^S",
    "B6", "B^R6^R", "B^R6^S", "B^S6^R", "B^S6^S",
    "B7", "B^R7^R", "B^R7^S", "B^S7^R", "B^S7^S",
    "B8", "B^R8", "B^S8",
    "B9", "B^R9^R", "B^R9^S", "B^S9^R", "B^S9^S",
    "B10", "B^R10^R", "B^R10^S", "B^S10^R", "B^S10^S",
    "B11",
    "B12", "B^R12^R", "B^R12^S", "B^S12^R", "B^S12^S",
    "B13", "B^R13", "B^S13",
]
PK_C_SUBTYPES = [
    "C1",
    "C2",
    "C4",
    "C7",
    "C13",
]
PK_D_SUBTYPES = [
    "D1",
    "D2", "D2^R", "D2^S",
    "D3",
    "D4", "D4^R", "D4^S",
    "D5", "D5^R", "D5^S",
    "D6", "D6^R", "D6^S",
    "D7", "D7^R", "D7^S",
    "D8", "D^R8", "D^S8",
    "D10", "D^R10^R", "D^R10^S", "D^S10^R", "D^S10^S",
    "D11",
    "D14",
    "D15",
    "D16", "D16^R", "D16^S",
    "D17", "D^R17", "D^S17",
]
HARDCODED_PK_SCORING = []
HARDCODED_PK_SCORING.extend([("A", pk_subtype, 1.0) for pk_subtype in PK_A_SUBTYPES])
HARDCODED_PK_SCORING.extend([("B", pk_subtype, 1.0) for pk_subtype in PK_B_SUBTYPES])
HARDCODED_PK_SCORING.extend([("C", pk_subtype, 1.0) for pk_subtype in PK_C_SUBTYPES])
HARDCODED_PK_SCORING.extend([("D", pk_subtype, 1.0) for pk_subtype in PK_D_SUBTYPES])


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db-path", type=str, required=True)
    parser.add_argument("--vocab-path", type=str, required=False)
    return parser.parse_args()


def main() -> None:
    args = cli()

    db_path: Path = Path(args.db_path).expanduser()

    # Setup fingerprinter
    rules = RuleSet.load_default()
    if not args.vocab_path:
        vocab_tokens = set()
        for rule in tqdm(rules.matching_rules, desc="building vocab"):
            rule_name = rule.name
            vocab_tokens.add(rule_name)
            rule_pseudonyms = rule.pseudonyms
            vocab_tokens.update(set(rule_pseudonyms))
        vocab = Vocabulary(vocab_tokens)
    else:
        vocab_path: Path = Path(args.vocab_path)
        vocab = Vocabulary.from_file(vocab_path)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)

    # Create scoring matrix dynamically, only match on diagonal
    records = []
    for rule in tqdm(rules.matching_rules, desc="extracting records"):
        rule_name = rule.name
        rule_smiles = rule.smiles
        records.append((rule_name, rule_smiles))
    scoring_matrix = create_tanimoto_scoring_matrix(
        records,
        radius=2,
        num_bits=2048,
        stereochemistry=False,
        self_score_tokens=[TOKEN_UNK, "A", "B", "C", "D"],
        self_score=1.0,
        hardcoded_scores=HARDCODED_PK_SCORING,
    )
    aligner = setup_aligner(scoring_matrix, mode="global")
    converter = Converter(
        to_identifier=lambda s: scoring_matrix.alphabet.index(s),
        from_identifier=lambda i: scoring_matrix.alphabet[i],
    )

    # Another test, very minimal fingerprint
    per_monomer_tokens = [
        ["PK_B", "PK"],
        ["PK_B", "PK"],
        ["PK_A", "PK"],
        ["PK_D", "PK"],
        ["PK_B", "PK"],
        ["PK_B", "PK"],
    ]
    fingerprint = fingerprinter.encode(per_monomer_tokens)
    with open_database(db_path, read_only=True) as db:
        results = db.closest(fingerprint, limit=1000)

        targets = [result.entry.primary_sequence for result in results]
        reranked = rerank(
            query=["B", "B", "A", "D", "B", "B"],
            targets=targets,
            aligner=aligner,
            converter=converter,
        )

        reranked_results = []
        for result, (score, inverted) in zip(results, reranked, strict=True):
            primary_sequence = result.entry.primary_sequence

            if inverted:
                primary_sequence = primary_sequence[::-1]

            reranked_results.append(
                (
                    score,
                    inverted,
                    result,
                    primary_sequence,
                )
            )

        reranked_results.sort(key=lambda x: x[0], reverse=True)

        top_n = 10
        for score, inverted, result, primary_sequence in reranked_results[:top_n]:
            print(score, result.similarity, result.entry.name)


if __name__ == "__main__":
    main()
