"""Shared helpers for the benchmark pipeline: discovery context, BGC<->compound
linking, and plotting style.

The fingerprint/alignment recipe in `build_discovery_context` must match
gui/src/server/routes/discovery.py's `_build_context` exactly (vocabulary, n_hashes=2,
Tanimoto scoring matrix radius/num_bits); see that module's docstring. Reproduced
here rather than imported, since discovery.py pulls in the Flask app (blueprints,
routes.database, routes.queue) this standalone pipeline has no use for.
"""

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from retromol.model.rules import MatchingRule, RuleSet
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import Converter
from retromol_alignment.scoring import HARDCODED_PK_SCORING, create_tanimoto_scoring_matrix
from retromol_database.duckdb import RetroMolDuckDB, FINGERPRINT_SIZE
from retromol_fingerprint.fingerprint import TOKEN_LINK, TOKEN_UNK, Fingerprinter, Vocabulary

PK_GROUP_TOKENS = ("PK_A", "PK_B", "PK_C", "PK_D")

# Both compound and bgc entries sourced from MIBiG carry a `mibig_url(accession,
# version)` URL (see database/scripts/common.py). The accession is the only
# queryable signal linking a BGC entry to its compound(s), since the database schema
# doesn't persist mibig_accession as its own column (see MIBIG_URL_TEMPLATE).
_MIBIG_ACCESSION_RE = re.compile(r"/repository/([A-Za-z0-9]+)\.")


def mibig_accession_from_url(url: str | None) -> str | None:
    if not url:
        return None
    m = _MIBIG_ACCESSION_RE.search(url)
    return m.group(1) if m else None


@dataclass(frozen=True)
class DiscoveryContext:
    ruleset: RuleSet
    name_to_rule: dict[str, MatchingRule]
    fingerprinter: Fingerprinter
    aligner: Any
    converter: Converter
    alphabet: frozenset[str]


def build_discovery_context(
    reaction_rules_path: str | Path | None = None,
    matching_rules_path: str | Path | None = None,
) -> DiscoveryContext:
    rules = RuleSet.load_from_files(
        reaction_rules_path=reaction_rules_path or None,
        matching_rules_path=matching_rules_path or None,
        match_stereochemistry=False,
    )

    vocab_tokens: set[str] = set()
    name_to_rule: dict[str, MatchingRule] = {}
    records: list[tuple[str, str]] = []

    for rule in rules.matching_rules:
        vocab_tokens.add(rule.name)
        vocab_tokens.update(rule.pseudonyms)
        name_to_rule.setdefault(rule.name, rule)
        records.append((rule.name, rule.smiles))

    vocab = Vocabulary(vocab_tokens)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)

    scoring_matrix = create_tanimoto_scoring_matrix(
        records,
        radius=2,
        num_bits=2048,
        stereochemistry=False,
        self_score_tokens=[TOKEN_UNK, TOKEN_LINK, *PK_GROUP_TOKENS],
        self_score=1.0,
        hardcoded_scores=HARDCODED_PK_SCORING,
    )
    aligner = setup_aligner(scoring_matrix, mode="global")
    alphabet = scoring_matrix.alphabet
    converter = Converter(
        to_identifier=lambda s: alphabet.index(s),
        from_identifier=lambda i: alphabet[i],
    )

    return DiscoveryContext(
        ruleset=rules,
        name_to_rule=name_to_rule,
        fingerprinter=fingerprinter,
        aligner=aligner,
        converter=converter,
        alphabet=frozenset(alphabet),
    )


def per_monomer_tokens(name: str, ctx: DiscoveryContext) -> list[str]:
    if name == TOKEN_LINK:
        return [TOKEN_LINK]
    if name in PK_GROUP_TOKENS:
        return [name, "PK"]
    rule = ctx.name_to_rule.get(name)
    if rule is None:
        return []
    tokens = {rule.name}
    tokens.update(rule.pseudonyms)
    return list(tokens)


def normalize_for_alignment(name: str, ctx: DiscoveryContext) -> str:
    if name == "X" or name not in ctx.alphabet:
        return TOKEN_UNK
    return name


def fingerprint_for_sequence(names: list[str], ctx: DiscoveryContext):
    tokens = [per_monomer_tokens(n, ctx) for n in names if n != TOKEN_LINK]
    return ctx.fingerprinter.encode(tokens)


def open_db(db_path: str | Path, read_only: bool = True) -> RetroMolDuckDB:
    return RetroMolDuckDB.open(db_path, read_only=read_only)


def bgc_biosynthetic_label_sets(con) -> dict[str, frozenset[str]]:
    """entry_id -> the full, exact set of biosynthetic_class_annotations labels for
    that bgc entry (e.g. {"Polyketide"}, or {"Polyketide", "Nonribosomal peptide"}
    for a PKS/NRPS hybrid). `con` is a duckdb connection (RetroMolDuckDB.con or a raw
    duckdb.connect(...) both work)."""
    labels_by_entry: dict[str, set[str]] = {}
    for entry_id, label in con.execute("SELECT entry_id, label FROM biosynthetic_class_annotations").fetchall():
        labels_by_entry.setdefault(entry_id, set()).add(label)
    return {entry_id: frozenset(labels) for entry_id, labels in labels_by_entry.items()}


def exact_label_group_accessions(
    con, links: dict[str, dict], group_label_sets: dict[str, frozenset[str]]
) -> dict[str, set[str]]:
    """{group_key: {accession, ...}} for benchmark/config.yaml's mibig_highlight_classes,
    where an accession belongs to a group only if its representative bgc entry's
    biosynthetic_class label set is EXACTLY that group's `biosynthetic_class_labels`
    (not a superset) -- so a BGC labeled e.g. Polyketide + Terpene lands in neither
    the "polyketide" nor any hybrid group, only in "all". The representative bgc
    entry is `links[accession]["bgc_entry_ids"][0]` (already sorted; same one
    retrieval_engine.build_eval_bgcs uses as that accession's query), so group
    membership always refers to the same bgc entry actually evaluated elsewhere.
    """
    label_sets = bgc_biosynthetic_label_sets(con)
    groups: dict[str, set[str]] = {key: set() for key in group_label_sets}

    for accession, ids in links.items():
        if not ids["bgc_entry_ids"]:
            continue
        bgc_labels = label_sets.get(ids["bgc_entry_ids"][0], frozenset())
        for key, required in group_label_sets.items():
            if bgc_labels == required:
                groups[key].add(accession)

    return groups


def figure_style() -> None:
    """Shared matplotlib backend/rcParams for every figureN_*.py script.

    Snakemake's `run:` rules execute in worker threads of the main process (not
    subprocesses) when multiple independent rules run under `--cores N`, so with
    figures 1-4 all runnable in parallel, matplotlib's default interactive backend
    (MacOSX on a Mac) crashes with "Cannot create a GUI FigureManager outside of the
    main thread" for every figure but whichever happened to run first. Forcing the
    non-interactive Agg backend (all these scripts only ever call savefig, never
    show()) fixes this regardless of which thread it runs on.
    """
    import matplotlib

    matplotlib.use("Agg", force=True)

    matplotlib.rcParams.update(
        {
            "figure.dpi": 150,
            "savefig.dpi": 150,
            "font.size": 10,
            "axes.titlesize": 11,
            "axes.labelsize": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "legend.frameon": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )
