"""Step 8: turn parsed BGC readouts into database entries.

Each antiSMASH region becomes one "bgc" entry. bgc_primary_sequence maps every
module's predicted substrate onto the same matching-rule vocabulary compound
primary sequences use (PKS modules resolve to a reduction-level pseudonym like
"PK_A", NRPS modules resolve by matching PARAS' predicted substrate structurally
against the ruleset) -- see retromol_antismash.modules for why the two module
types resolve differently. That shared vocabulary is what makes a BGC entry's
fingerprint and primary sequence comparable to a compound's.

The MIBiG URL needs an accession's version, which parse_gbks.py can no longer
provide (MIBiG 4.0's GBKs dropped the ACCESSION.VERSION suffix) -- so it's
looked up here instead, from the same JSON-derived mibig_versions.json
extract_mibig_compounds.py produces for load_compounds.py.

BGCs are deduplicated per source .gbk file via `file_hash` (parse_gbks.py's
sha256 of the whole file's raw text): a region's entry id is derived from that
hash plus its own readout id, and a file whose hash is already stored (checked
via RetroMolDuckDB.bgc_content_hash_exists) is skipped outright, so rerunning
the pipeline over an already-ingested file doesn't redo any fingerprinting work.
"""

import argparse
import hashlib
import json
import logging
from pathlib import Path

from tqdm import tqdm

from common import build_fingerprint_context, load_ruleset, mibig_url, phylogeny_from_organism_name
from retromol_antismash.modules import LinearReadout, bgc_primary_sequence
from retromol_database.duckdb import RetroMolDuckDB
from taxonomy import TaxonomyDB, resolve_phylogeny

log = logging.getLogger(__name__)


def run(
    readouts_path: str | Path,
    db_path: str | Path,
    reaction_rules_path: str | Path | None,
    matching_rules_path: str | Path | None,
    mibig_versions_path: str | Path,
    mibig_annotations_path: str | Path,
    match_stereochemistry: bool = False,
    include_raw_gbk: bool = True,
    taxdump_dir: str | Path | None = None,
) -> None:
    ruleset = load_ruleset(reaction_rules_path, matching_rules_path, match_stereochemistry)
    _, fingerprinter = build_fingerprint_context(ruleset)

    taxdb = TaxonomyDB.load(taxdump_dir) if taxdump_dir else None

    with open(mibig_versions_path) as fh:
        versions: dict[str, str] = json.load(fh)

    with open(mibig_annotations_path) as fh:
        annotations: dict[str, dict] = json.load(fh)

    added = 0
    skipped = 0
    skipped_existing_file = 0

    db = RetroMolDuckDB.open(db_path)
    try:
        with open(readouts_path) as fh:
            pbar = tqdm(fh, desc="load_bgcs", unit="region")
            for line in pbar:
                line = line.strip()
                if not line:
                    continue

                entry = json.loads(line)
                file_hash = entry.get("file_hash")

                if file_hash and db.bgc_content_hash_exists(file_hash):
                    skipped_existing_file += 1
                    pbar.set_postfix(added=added, skipped=skipped, skipped_existing=skipped_existing_file)
                    continue

                readout = LinearReadout.from_dict(entry["readout"])
                names, tokens = bgc_primary_sequence(readout, ruleset)

                if not names:
                    skipped += 1
                    pbar.set_postfix(added=added, skipped=skipped, skipped_existing=skipped_existing_file)
                    continue

                fp = fingerprinter.encode(tokens)
                accession = entry.get("accession")
                name = f"{accession} ({readout.id})" if accession else readout.id
                url = mibig_url(accession, versions.get(accession)) if accession else None
                entry_id = hashlib.sha256(f"{file_hash}:{readout.id}".encode("utf-8")).hexdigest()

                db.add_entry(
                    entry_id=entry_id,
                    name=name,
                    database_name="MIBiG",
                    url=url,
                    raw=entry.get("raw_gbk") if include_raw_gbk else None,
                    entry_type="bgc",
                    primary_sequence=names,
                    fingerprint=fp,
                    content_hash=file_hash,
                )

                record = annotations.get(accession) if accession else None
                if record:
                    fallback_type, fallback_genus, fallback_species = phylogeny_from_organism_name(
                        record.get("organism_name")
                    )
                    resolution = resolve_phylogeny(
                        taxdb,
                        ncbi_tax_id=record.get("ncbi_tax_id"),
                        genus=fallback_genus,
                        species=fallback_species,
                        fallback_type_label=fallback_type,
                    )
                    db.add_phylogeny_annotation(
                        entry_id,
                        type_label=resolution.type_label,
                        type_taxid=resolution.type_taxid,
                        genus=resolution.genus,
                        genus_taxid=resolution.genus_taxid,
                        species=resolution.species,
                        species_taxid=resolution.species_taxid,
                    )
                    # biosynthetic_class describes this BGC's own biosynthesis machinery
                    # (PKS/NRPS/RiPP/...) -- a gene-cluster property, not populated for
                    # compounds (see RetroMolDuckDB.add_biosynthetic_class_annotation /
                    # load_compounds.py).
                    for chemical_class in record.get("biosyn_class") or []:
                        if chemical_class:
                            db.add_biosynthetic_class_annotation(entry_id, str(chemical_class))

                added += 1
                pbar.set_postfix(added=added, skipped=skipped, skipped_existing=skipped_existing_file)
    finally:
        db.close()

    log.info(
        "load_bgcs: added=%d skipped=%d skipped_existing_file=%d",
        added, skipped, skipped_existing_file,
    )


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--readouts", required=True)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--rxn-rules", default=None)
    ap.add_argument("--mxn-rules", default=None)
    ap.add_argument("--mibig-versions", required=True)
    ap.add_argument("--mibig-annotations", required=True)
    ap.add_argument("--match-stereochemistry", action="store_true")
    ap.add_argument("--no-raw-gbk", action="store_true")
    ap.add_argument("--taxdump-dir", default=None, help="dir with NCBI names.dmp/nodes.dmp (skips taxid resolution if omitted)")
    args = ap.parse_args()

    run(
        readouts_path=args.readouts,
        db_path=args.db_path,
        reaction_rules_path=args.rxn_rules,
        matching_rules_path=args.mxn_rules,
        mibig_versions_path=args.mibig_versions,
        mibig_annotations_path=args.mibig_annotations,
        match_stereochemistry=args.match_stereochemistry,
        include_raw_gbk=not args.no_raw_gbk,
        taxdump_dir=args.taxdump_dir,
    )


if __name__ == "__main__":
    main()
