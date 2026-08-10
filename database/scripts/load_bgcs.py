"""Step 8: turn parsed BGC readouts into database entries.

Each antiSMASH region becomes one "bgc" entry. bgc_primary_sequence maps every
module's predicted substrate onto the same matching-rule vocabulary compound
primary sequences use (PKS modules resolve to a reduction-level pseudonym like
"PK_A", NRPS modules resolve by matching PARAS' predicted substrate structurally
against the ruleset) -- see retromol_antismash.modules for why the two module
types resolve differently. That shared vocabulary is what makes a BGC entry's
fingerprint and primary sequence comparable to a compound's.
"""

import argparse
import json
import logging
from pathlib import Path

from common import build_fingerprint_context, load_ruleset
from retromol_antismash.modules import LinearReadout, bgc_primary_sequence
from retromol_database.duckdb import RetroMolDuckDB

log = logging.getLogger(__name__)


def run(
    readouts_path: str | Path,
    db_path: str | Path,
    reaction_rules_path: str | Path | None,
    matching_rules_path: str | Path | None,
    match_stereochemistry: bool = False,
    include_raw_gbk: bool = True,
) -> None:
    ruleset = load_ruleset(reaction_rules_path, matching_rules_path, match_stereochemistry)
    _, fingerprinter = build_fingerprint_context(ruleset)

    added = 0
    skipped = 0

    db = RetroMolDuckDB.open(db_path)
    try:
        with open(readouts_path) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue

                entry = json.loads(line)
                readout = LinearReadout.from_dict(entry["readout"])
                names, tokens = bgc_primary_sequence(readout, ruleset)

                if not names:
                    skipped += 1
                    continue

                fp = fingerprinter.encode(tokens)
                name = f"{entry['accession']} ({readout.id})" if entry.get("accession") else readout.id

                db.add_entry(
                    name=name,
                    url=entry.get("url"),
                    raw=entry.get("raw_gbk") if include_raw_gbk else None,
                    entry_type="bgc",
                    primary_sequence=names,
                    fingerprint=fp,
                )
                added += 1
    finally:
        db.close()

    log.info("load_bgcs: added=%d skipped=%d", added, skipped)


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--readouts", required=True)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--rxn-rules", default=None)
    ap.add_argument("--mxn-rules", default=None)
    ap.add_argument("--match-stereochemistry", action="store_true")
    ap.add_argument("--no-raw-gbk", action="store_true")
    args = ap.parse_args()

    run(
        readouts_path=args.readouts,
        db_path=args.db_path,
        reaction_rules_path=args.rxn_rules,
        matching_rules_path=args.mxn_rules,
        match_stereochemistry=args.match_stereochemistry,
        include_raw_gbk=not args.no_raw_gbk,
    )


if __name__ == "__main__":
    main()
