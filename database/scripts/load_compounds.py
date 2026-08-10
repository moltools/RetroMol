"""Steps 4 & 6: turn parsed compound results into database entries.

For each RetroMol Result, every reconstructed candidate primary sequence
(retromol_synthesis.reconstruction.reconstruct_linear_readout) becomes its own
"compound" entry -- ambiguous parses intentionally produce multiple queryable
entries, the same convention the webapp uses for an uploaded compound with more
than one candidate reading.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Literal

from common import build_fingerprint_context, find_key_ci, load_ruleset, mibig_url, npatlas_url, per_monomer_tokens
from retromol.model.result import Result
from retromol_database.duckdb import RetroMolDuckDB
from retromol_synthesis.reconstruction import reconstruct_linear_readout

log = logging.getLogger(__name__)


def _npatlas_name_and_url(props: dict) -> tuple[str | None, str | None]:
    npaid_key = find_key_ci(props, ["npaid"])
    npaid = props.get(npaid_key) if npaid_key else None

    name_key = find_key_ci(props, ["original_name", "compound_name", "name"])
    name = props.get(name_key) if name_key else None

    return name, npatlas_url(npaid)


def _mibig_name_and_url(props: dict, versions: dict[str, str]) -> tuple[str | None, str | None]:
    accession = props.get("mibig_accession")
    version = versions.get(accession) if accession else None
    name = props.get("name")
    return name, mibig_url(accession, version)


def run(
    results_path: str | Path,
    db_path: str | Path,
    source: Literal["npatlas", "mibig"],
    reaction_rules_path: str | Path | None,
    matching_rules_path: str | Path | None,
    match_stereochemistry: bool = False,
    mibig_versions_path: str | Path | None = None,
) -> None:
    ruleset = load_ruleset(reaction_rules_path, matching_rules_path, match_stereochemistry)
    name_to_rule, fingerprinter = build_fingerprint_context(ruleset)

    versions: dict[str, str] = {}
    if source == "mibig" and mibig_versions_path is not None:
        with open(mibig_versions_path) as fh:
            versions = json.load(fh)

    added = 0
    skipped = 0

    db = RetroMolDuckDB.open(db_path)
    try:
        with open(results_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue

                result = Result.from_dict(json.loads(line))
                props = result.submission.props or {}

                if source == "npatlas":
                    name, url = _npatlas_name_and_url(props)
                else:
                    name, url = _mibig_name_and_url(props, versions)

                name = name or result.submission.name or result.submission.inchikey

                for reconstruction in reconstruct_linear_readout(result):
                    names = [n for n, _tags in reconstruction.primary_sequence]
                    if not names:
                        skipped += 1
                        continue

                    tokens = [per_monomer_tokens(n, name_to_rule) for n in names]
                    fp = fingerprinter.encode(tokens)

                    db.add_entry(
                        name=name,
                        url=url,
                        raw=result.submission.smiles,
                        entry_type="compound",
                        primary_sequence=names,
                        fingerprint=fp,
                    )
                    added += 1
    finally:
        db.close()

    log.info("load_compounds[%s]: added=%d skipped=%d", source, added, skipped)


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--results", required=True)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--source", choices=["npatlas", "mibig"], required=True)
    ap.add_argument("--rxn-rules", default=None)
    ap.add_argument("--mxn-rules", default=None)
    ap.add_argument("--match-stereochemistry", action="store_true")
    ap.add_argument("--mibig-versions", default=None, help="required when --source=mibig")
    args = ap.parse_args()

    run(
        results_path=args.results,
        db_path=args.db_path,
        source=args.source,
        reaction_rules_path=args.rxn_rules,
        matching_rules_path=args.mxn_rules,
        match_stereochemistry=args.match_stereochemistry,
        mibig_versions_path=args.mibig_versions,
    )


if __name__ == "__main__":
    main()
