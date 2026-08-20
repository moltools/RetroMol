"""Shared helpers for the database-construction pipeline.

The fingerprint recipe here (vocabulary = every matching rule's name + pseudonyms,
Fingerprinter with n_bits=FINGERPRINT_SIZE, n_hashes=2) must match
gui/src/server/routes/discovery.py's `_build_context` exactly -- that's what the
webapp uses to encode a query at search time. Deviating here would silently make
every fingerprint stored by this pipeline incomparable to a live query.
"""

from multiprocessing import Pool
from pathlib import Path
from typing import Any, Iterable, Iterator

from rdkit import RDLogger

from retromol.io.streaming import ResultEvent, _init_worker, _process_compound, _task_buffered_iterator
from retromol.model.result import Result
from retromol.model.rules import MatchingRule, RuleSet
from retromol_database.duckdb import FINGERPRINT_SIZE
from retromol_fingerprint.fingerprint import TOKEN_LINK, Fingerprinter, Vocabulary

# Silences RDKit's kekulization/valence/etc. warnings in *this* (single) process --
# every pipeline script imports common, so this alone covers create_db.py,
# load_compounds.py, load_bgcs.py, and parse_gbks.py's main process. It does NOT
# reliably reach parse_compounds.py's multiprocessing workers: those are spawned
# fresh by retromol.io.streaming.run_retromol_stream's own Pool, with its own
# initializer, and whether a fresh worker re-runs this module-level call at all
# depends on whether the *original* entry point was this script directly (works)
# or something else re-importing it, like Snakemake's generated run: script
# (doesn't) -- see run_retromol_stream_quiet below for the reliable fix.
RDLogger.DisableLog("rdApp.*")

# Group-level PKS pseudo-tokens a BGC's PKS module resolves to (see
# retromol_antismash.modules.PKSExtenderUnit / module_primary_sequence_tokens).
# Not matching-rule names themselves, but already part of the fingerprint
# vocabulary as pseudonyms of every rule at that reduction level.
PK_GROUP_TOKENS = ("PK_A", "PK_B", "PK_C", "PK_D")

MIBIG_URL_TEMPLATE = "https://mibig.secondarymetabolites.org/repository/{accession}.{version}/index.html#r1c1"
NPATLAS_URL_TEMPLATE = "https://www.npatlas.org/explore/compounds/{npaid}"


def load_ruleset(
    reaction_rules_path: str | Path | None,
    matching_rules_path: str | Path | None,
    match_stereochemistry: bool = False,
) -> RuleSet:
    """Load a RuleSet, falling back to RetroMol's bundled default rules when a path is None/empty."""
    return RuleSet.load_from_files(
        reaction_rules_path=reaction_rules_path or None,
        matching_rules_path=matching_rules_path or None,
        match_stereochemistry=match_stereochemistry,
    )


def build_fingerprint_context(ruleset: RuleSet) -> tuple[dict[str, MatchingRule], Fingerprinter]:
    """Build the (name_to_rule, fingerprinter) pair used to fingerprint primary sequences."""
    vocab_tokens: set[str] = set()
    name_to_rule: dict[str, MatchingRule] = {}

    for rule in ruleset.matching_rules:
        vocab_tokens.add(rule.name)
        vocab_tokens.update(rule.pseudonyms)
        name_to_rule.setdefault(rule.name, rule)

    vocab = Vocabulary(vocab_tokens)
    fingerprinter = Fingerprinter(vocab, n_bits=FINGERPRINT_SIZE, n_hashes=2)
    return name_to_rule, fingerprinter


def per_monomer_tokens(name: str, name_to_rule: dict[str, MatchingRule]) -> list[str]:
    """Fingerprinting token list for one compound primary sequence block (mirrors discovery.py's `_per_monomer_tokens`).

    TOKEN_LINK is not a building block -- it just joins two merged paths -- and is
    filtered out by callers before fingerprinting (see load_compounds.py), so it's
    never actually looked up here. Handled explicitly anyway so a stray call can't
    silently fall through to the empty-tokens/TOKEN_UNK path instead.
    """
    if name == TOKEN_LINK:
        return [TOKEN_LINK]

    if name in PK_GROUP_TOKENS:
        return [name, "PK"]

    rule = name_to_rule.get(name)
    if rule is None:
        return []

    tokens = {rule.name}
    tokens.update(rule.pseudonyms)
    return list(tokens)


def find_key_ci(record: dict[str, Any], substrings: list[str]) -> str | None:
    """Find the first key in `record` whose lowercased form contains any of `substrings`."""
    for key in record:
        lowered = str(key).lower()
        if any(sub in lowered for sub in substrings):
            return key
    return None


def mibig_url(accession: str | None, version: str | None) -> str | None:
    if not accession or not version:
        return None
    return MIBIG_URL_TEMPLATE.format(accession=accession, version=version)


def npatlas_url(npaid: str | None) -> str | None:
    if not npaid:
        return None
    return NPATLAS_URL_TEMPLATE.format(npaid=npaid)


def primary_sequence_from_result(result: Result) -> list[str]:
    """
    The single primary sequence for a parsed compound, read directly off
    `result.linear_readout` -- no backbone reconstruction involved, that's a
    display-only concern this pipeline has no use for. `result.linear_readout` is
    already computed by retromol.pipelines.parsing.run_retromol (it's just a field
    on Result); every path through the molecule -- including single-node paths for
    tailoring events that don't connect to any chain, e.g. glycosylation/methylation
    (AssemblyGraph only keeps C-C/C-N bonds as "connections", so a sugar attached via
    a glycosidic C-O-C linkage always ends up disconnected from the main chain) --
    is merged into one sequence (longest first, ties broken lexicographically, joined
    by TOKEN_LINK -- see `LinearReadout.primary_sequence`), and that single sequence
    becomes this compound's one db entry (see load_compounds.py). Nothing found in
    the assembly graph is dropped. An unidentified node is named "X", the same
    convention used everywhere else in RetroMol.

    :param result: a parsed RetroMol Result
    :return: the single primary sequence
    """
    return result.linear_readout.primary_sequence()


def _init_worker_quiet(ruleset: RuleSet) -> None:
    """Worker-process initializer: set up the ruleset global exactly like retromol.io.streaming's own
    _init_worker does, then also disable RDKit logging -- unlike a module-level RDLogger call, this is
    guaranteed to run once per worker process no matter how the pool was launched."""
    _init_worker(ruleset)
    RDLogger.DisableLog("rdApp.*")


def run_retromol_stream_quiet(
    ruleset: RuleSet,
    row_iter: Iterable[dict[str, Any]],
    smiles_col: str = "smiles",
    workers: int = 1,
    batch_size: int = 2000,
    pool_chunksize: int = 50,
    maxtasksperchild: int = 2000,
) -> Iterator[ResultEvent]:
    """
    Drop-in replacement for retromol.io.streaming.run_retromol_stream that also
    disables RDKit's C-level logging inside every worker process. Reuses that
    module's own batching/worker-task functions -- only the Pool's initializer
    differs (see _init_worker_quiet).
    """
    with Pool(
        processes=workers,
        initializer=_init_worker_quiet,
        initargs=(ruleset,),
        maxtasksperchild=maxtasksperchild,
    ) as pool:
        for task_batch in _task_buffered_iterator(row_iter, smiles_col=smiles_col, batch_size=batch_size):
            for serialized, err in pool.imap_unordered(_process_compound, task_batch, chunksize=pool_chunksize):
                yield ResultEvent(serialized, err)


def split_accession_version(record_id: str) -> tuple[str, str | None]:
    """
    Split a GenBank-style "ACCESSION.VERSION" id (e.g. "BGC0000001.5") in two.

    :param record_id: the record id, as set on retromol_antismash Region/LinearReadout.id
    :return: (accession, version) -- version is None if record_id has no "." suffix
    """
    if "." in record_id:
        accession, version = record_id.rsplit(".", 1)
        return accession, version
    return record_id, None
