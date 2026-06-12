import argparse
from pathlib import Path

from tqdm import tqdm

from retromol_database.duckdb import FINGERPRINT_SIZE, open_database
from retromol.model.rules import RuleSet
from retromol_fingerprint.fingerprint import TOKEN_UNK, Vocabulary, Fingerprinter
from retromol_alignment.scoring import create_tanimoto_scoring_matrix, HARDCODED_PK_SCORING
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import Converter
from retromol_alignment.ranking import rerank


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
    # per_monomer_tokens = [
    #     ["PK_A", "PK"],
    #     ["cysteine"],
    #     ["PK_C", "PK"],
    #     ["PK_B", "PK"],
    #     ["PK_B", "PK"],
    #     ["PK_D", "PK"],
    #     ["PK_D", "PK"],
    #     ["PK_B", "PK"],
    #     ["PK_A", "PK"],
    #     ["PK_D", "PK"],
    # ]
    # per_monomer_tokens = [
    #     ["tryptophan"],
    #     ["asparagine"],
    #     ["aspartic acid"],
    #     ["threonine"],
    #     ["glycine"],
    #     ["ornithine"],
    #     ["aspartic acid"],
    #     ["alanine"],
    #     ["aspartic acid"],
    #     ["glycine"],
    #     ["serine"],
    #     ["3-methylglutamic acid"],
    #     ["kynurenine"]
    # ]
    # per_monomer_tokens = [
    #     # ["acetic acid"],
    #     # ["PK_D", "D1"],
    #     # ["PK_D", "D1"],
    #     # ["PK_D", "D1"],
    #     # ["PK_D", "D1"],
    #     # ["PK_D", "D1"],
    #     # ["PK_C", "C1"],
    #     ["N-(5-aminopentyl)hydroxylamine", "siderophore"],
    #     ["butanedioic acid", "siderophore"],
    #     ["N-(5-aminopentyl)hydroxylamine", "siderophore"],
    #     ["serine"],
    #     ["salicylic acid"],
    # ]
    fingerprint = fingerprinter.encode(per_monomer_tokens)
    with open_database(db_path, read_only=True) as db:
        results = db.closest(fingerprint, limit=1000)

        targets = [result.entry.primary_sequence for result in results]
        reranked = rerank(
            query=["B", "B", "A", "D", "B", "B"],
            # query=["A", "cysteine", "C", "B", "B", "D", "D", "B", "A", "D"],
            # query=[
            #     "tryptophan",
            #     "asparagine",
            #     "aspartic acid",
            #     "threonine",
            #     "glycine",
            #     "ornithine",
            #     "aspartic acid",
            #     "alanine",
            #     "aspartic acid",
            #     "glycine",
            #     "serine",
            #     "3-methylglutamic acid",
            #     "kynurenine",
            # ],
            # query=[
            #     # "acetic acid",
            #     # "D",
            #     # "D",
            #     # "D",
            #     # "D",
            #     # "D",
            #     # "C",
            #     "N-(5-aminopentyl)hydroxylamine",
            #     "butanedioic acid",
            #     "N-(5-aminopentyl)hydroxylamine",
            #     "serine",
            #     "salicylic acid",
            # ],
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

        top_n = 20
        for score, inverted, result, primary_sequence in reranked_results[:top_n]:
            print(score, result.similarity, result.entry.name, result.entry.raw)


if __name__ == "__main__":
    main()
