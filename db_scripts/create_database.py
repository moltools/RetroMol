import argparse
from pathlib import Path

from tqdm import tqdm

from retromol.io.json import iter_json
from retromol_database.duckdb import FINGERPRINT_SIZE, create_database, open_database
from retromol.model.result import Result
from retromol.model.rules import RuleSet
from retromol.model.reaction_graph import MolNode
from retromol_fingerprint.fingerprint import TOKEN_UNK, Vocabulary, Fingerprinter
from retromol_alignment.scoring import create_tanimoto_scoring_matrix
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import Converter
from retromol_alignment.ranking import rerank


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db-path", type=str, required=True)
    parser.add_argument("--npatlas-retromol-path", type=str, required=True)
    parser.add_argument("--create-new", action="store_true")
    return parser.parse_args()


def parse_molnode_path(node_path: list[MolNode]) -> tuple[list[str], list[list[str]]]:
    primary_sequence = []
    per_monomer_tokens = []

    for node in node_path:
        if node.identified and node.identity:
            identity = node.identity.matched_rule
            monomer_tokens = set()
            monomer_tokens.add(identity.name)
            monomer_tokens.update(set(identity.pseudonyms))
            primary_sequence.append(identity.name)
            per_monomer_tokens.append(list(monomer_tokens))
        else:
            primary_sequence.append(TOKEN_UNK)
            per_monomer_tokens.append([])  # gets translated to unknown token inside encoding function

    return primary_sequence, per_monomer_tokens


def main() -> None:
    args = cli()

    db_path: Path = Path(args.db_path).expanduser()
    npatlas_retromol_path: str = args.npatlas_retromol_path

    # Setup fingerprinter
    rules = RuleSet.load_default()
    vocab_tokens = set()
    for rule in tqdm(rules.matching_rules, desc="building vocab"):
        rule_name = rule.name
        vocab_tokens.add(rule_name)
        rule_pseudonyms = rule.pseudonyms
        vocab_tokens.update(set(rule_pseudonyms))
    vocab = Vocabulary(vocab_tokens)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)

    # Add NPAtlas entries to database, if overwriting
    if args.create_new:
        with create_database(db_path, overwrite=args.create_new) as db:
            for d in tqdm(iter_json(npatlas_retromol_path, jsonl=True), desc="adding npatlas to db"):
                result = Result.from_dict(d)

                props = result.submission.props or {}
                npaid = props.get("npaid", "")
                url = f"https://www.npatlas.org/explore/compounds/{npaid}" if npaid else None
                name = props.get("compound_name", "")
                raw = props.get("smiles", "")

                for node_path in result.linear_readout.paths:
                    primary_sequence, per_monomer_tokens = parse_molnode_path(node_path)
                    fingerprint = fingerprinter.encode(per_monomer_tokens)

                    db.add_entry(
                        name=name,
                        url=url,
                        raw=raw,
                        entry_type="compound",
                        primary_sequence=primary_sequence,
                        fingerprint=fingerprint,
                    )

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
        self_score_tokens=[TOKEN_UNK],
        self_score=1.0,
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
            query=["B1", "B1", "A1", "D1", "B1", "B1"],
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
