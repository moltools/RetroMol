"""Build the MIBiG accession -> {bgc_entry_ids, compound_entry_ids} mapping.

The compiled database has no direct BGC<->compound foreign key (see benchmark/
config.yaml's mibig_highlight_classes comment and common.mibig_accession_from_url):
both entry types carry the same `mibig_url(accession, version)` URL as one of their
`entry_sources` rows, which is the only queryable link between a BGC and its
product compound(s). Re-derived here once and cached, rather than re-parsed by every
downstream script.

Also drops any bgc entry with fewer than `min_modules` primary-sequence blocks (see
benchmark/config.yaml's min_bgc_modules) -- filtering here, rather than per-figure,
keeps every figure's MIBiG-side sample consistent, since they all read this one
links.json.
"""

import argparse
import json
import logging
from collections import defaultdict
from pathlib import Path

import duckdb

from common import mibig_accession_from_url
from retromol_fingerprint.fingerprint import TOKEN_LINK

log = logging.getLogger(__name__)


def run(db_path: str | Path, output_path: str | Path, min_modules: int = 1) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    con = duckdb.connect(str(db_path), read_only=True)
    try:
        rows = con.execute(
            """
            SELECT e.type, s.entry_id, s.url
            FROM entry_sources s
            JOIN entries e ON e.id = s.entry_id
            WHERE s.database_name = 'MIBiG' AND s.url IS NOT NULL
            """
        ).fetchall()

        # Module count per bgc entry: primary_sequence length excluding TOKEN_LINK
        # (a join marker between merged paths, not a biosynthetic module).
        bgc_module_counts: dict[str, int] = {}
        if min_modules > 1:
            for entry_id, primary_sequence in con.execute(
                "SELECT id, primary_sequence FROM entries WHERE type = 'bgc'"
            ).fetchall():
                bgc_module_counts[entry_id] = sum(1 for tok in primary_sequence if tok != TOKEN_LINK)
    finally:
        con.close()

    links: dict[str, dict[str, set]] = defaultdict(lambda: {"bgc_entry_ids": set(), "compound_entry_ids": set()})
    n_bgc_dropped = 0

    for entry_type, entry_id, url in rows:
        accession = mibig_accession_from_url(url)
        if not accession:
            continue
        if entry_type == "bgc":
            if min_modules > 1 and bgc_module_counts.get(entry_id, 0) < min_modules:
                n_bgc_dropped += 1
                continue
            links[accession]["bgc_entry_ids"].add(entry_id)
        else:
            links[accession]["compound_entry_ids"].add(entry_id)

    serializable = {
        accession: {
            "bgc_entry_ids": sorted(ids["bgc_entry_ids"]),
            "compound_entry_ids": sorted(ids["compound_entry_ids"]),
        }
        for accession, ids in links.items()
    }

    with open(output_path, "w") as fh:
        json.dump(serializable, fh, indent=2, sort_keys=True)

    n_with_both = sum(
        1 for v in serializable.values() if v["bgc_entry_ids"] and v["compound_entry_ids"]
    )
    log.info(
        "link_bgc_compounds: %d MIBiG accessions, %d with both a bgc and >=1 compound entry "
        "(min_modules=%d dropped %d bgc entries) -> %s",
        len(serializable), n_with_both, min_modules, n_bgc_dropped, output_path,
    )


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-modules", type=int, default=1)
    args = ap.parse_args()
    run(db_path=args.db_path, output_path=args.output, min_modules=args.min_modules)


if __name__ == "__main__":
    main()
