"""Shared helpers for the database-construction pipeline.

The fingerprint recipe here (vocabulary = every matching rule's name + pseudonyms,
Fingerprinter with n_bits=FINGERPRINT_SIZE, n_hashes=2) must match
gui/src/server/routes/discovery.py's `_build_context` exactly -- that's what the
webapp uses to encode a query at search time. Deviating here would silently make
every fingerprint stored by this pipeline incomparable to a live query.
"""

from pathlib import Path
from typing import Any

from retromol.model.rules import MatchingRule, RuleSet
from retromol_database.duckdb import FINGERPRINT_SIZE
from retromol_fingerprint.fingerprint import Fingerprinter, Vocabulary

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
    """Fingerprinting token list for one compound primary-sequence block (mirrors discovery.py's `_per_monomer_tokens`)."""
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
