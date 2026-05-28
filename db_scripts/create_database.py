import argparse
from pathlib import Path

from tqdm import tqdm

from retromol.io.json import iter_json
from retromol_database.duckdb import FINGERPRINT_SIZE, create_database
from retromol.model.result import Result
from retromol.model.rules import RuleSet
from retromol.model.reaction_graph import MolNode
from retromol_fingerprint.fingerprint import TOKEN_UNK, Vocabulary, Fingerprinter


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db-path", type=str, required=True)
    parser.add_argument("--vocab-path", type=str, required=True)
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
    vocab_path: Path = Path(args.vocab_path).expanduser()
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
    vocab.write(vocab_path)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)

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


if __name__ == "__main__":
    main()
