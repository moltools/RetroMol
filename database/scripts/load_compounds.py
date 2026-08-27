"""Steps 4 & 6: turn parsed compound results into database entries.

For each RetroMol Result, every path found in result.linear_readout -- including a
compound's tailoring events (glycosylation, methylation -- anything that shows up as
its own disconnected single-node path) -- is merged into the one primary sequence
stored for that compound: longest path first, ties broken lexicographically, joined
by TOKEN_LINK (see common.primary_sequence_from_result). `raw` is the original input
SMILES.

Compounds are deduplicated on `result.submission.inchikey` (stereo-aware, computed by
retromol.model.submission.Submission with keep_stereo=True). The same molecule
appearing in both NPAtlas and MIBiG lands as one database entry with two source
records (see RetroMolDuckDB.add_entry) rather than two separate rows.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Literal

from tqdm import tqdm

from common import (
    build_fingerprint_context,
    clean_genus,
    clean_species_epithet,
    find_key_ci,
    load_ruleset,
    mibig_url,
    npatlas_url,
    per_monomer_tokens,
    phylogeny_from_organism_name,
    primary_sequence_from_result,
)
from retromol.model.result import Result
from retromol_database.duckdb import RetroMolDuckDB
from retromol_fingerprint.fingerprint import TOKEN_LINK
from taxonomy import TaxonomyDB, resolve_phylogeny

log = logging.getLogger(__name__)

DATABASE_NAMES = {"npatlas": "NPAtlas", "mibig": "MIBiG"}


def _npatlas_name_and_url(props: dict) -> tuple[str | None, str | None]:
    npaid_key = find_key_ci(props, ["npaid"])
    npaid = props.get(npaid_key) if npaid_key else None

    name_key = find_key_ci(props, ["original_name", "compound_name", "name"])
    name = props.get(name_key) if name_key else None

    return name, npatlas_url(npaid)


def _npatlas_phylogeny(props: dict) -> tuple[str | None, str | None, str | None]:
    """NPAtlas's SDF carries type/genus/species directly (unlike MIBiG's free-text
    organism_name) -- see the `origin_type`/`genus`/`origin_species` SDF properties.

    Both genus and origin_species are used as-is by NPAtlas's own curators, which
    turns out to include the same non-taxonomic placeholders MIBiG's free text does
    (bare "sp.", a genus name duplicated into the species field e.g. "Streptomyces
    sp.", "unidentified") -- cleaned the same way as MIBiG's organism_name parsing
    (see common.clean_genus/clean_species_epithet), so neither source's placeholder
    values end up stored as if they were real species-level identifications.
    """
    type_key = find_key_ci(props, ["origin_type"])
    type_label = props.get(type_key) if type_key else None

    genus_key = find_key_ci(props, ["genus"])
    genus = clean_genus(props.get(genus_key) if genus_key else None)

    species_key = find_key_ci(props, ["origin_species"])
    species_raw = props.get(species_key) if species_key else None
    species = clean_species_epithet(genus, species_raw)

    return (type_label or None), genus, species


def _mibig_name_and_url(props: dict, versions: dict[str, str]) -> tuple[str | None, str | None]:
    accession = props.get("mibig_accession")
    version = versions.get(accession) if accession else None
    name = props.get("name")
    return name, mibig_url(accession, version)


def _apply_mibig_annotations(
    db: RetroMolDuckDB, entry_id: str, props: dict, annotations: dict[str, dict], taxdb: TaxonomyDB | None
) -> None:
    accession = props.get("mibig_accession")
    record = annotations.get(accession) if accession else None
    if not record:
        return

    fallback_type, fallback_genus, fallback_species = phylogeny_from_organism_name(record.get("organism_name"))
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
    # biosynthetic_class describes the BGC's own biosynthesis machinery (PKS/NRPS/...),
    # not the compound structure -- populated in load_bgcs.py instead, not here.


def run(
    results_path: str | Path,
    db_path: str | Path,
    source: Literal["npatlas", "mibig"],
    reaction_rules_path: str | Path | None,
    matching_rules_path: str | Path | None,
    match_stereochemistry: bool = False,
    mibig_versions_path: str | Path | None = None,
    mibig_annotations_path: str | Path | None = None,
    taxdump_dir: str | Path | None = None,
    log_every: int = 1000,
) -> None:
    ruleset = load_ruleset(reaction_rules_path, matching_rules_path, match_stereochemistry)
    name_to_rule, fingerprinter = build_fingerprint_context(ruleset)

    taxdb = TaxonomyDB.load(taxdump_dir) if taxdump_dir else None

    versions: dict[str, str] = {}
    if source == "mibig" and mibig_versions_path is not None:
        with open(mibig_versions_path) as fh:
            versions = json.load(fh)

    annotations: dict[str, dict] = {}
    if source == "mibig" and mibig_annotations_path is not None:
        with open(mibig_annotations_path) as fh:
            annotations = json.load(fh)

    compounds = 0
    added = 0
    skipped = 0

    db = RetroMolDuckDB.open(db_path)
    try:
        with open(results_path) as fh:
            with tqdm(desc=f"load_compounds[{source}]", unit="cmpd") as pbar:
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

                    names = primary_sequence_from_result(result)

                    if not names:
                        skipped += 1
                    else:
                        # TOKEN_LINK only marks where two merged paths join -- it isn't
                        # a building block, so it's excluded from the fingerprint (an
                        # empty token list would otherwise silently add TOKEN_UNK mass).
                        tokens = [per_monomer_tokens(n, name_to_rule) for n in names if n != TOKEN_LINK]
                        fp = fingerprinter.encode(tokens)

                        db.add_entry(
                            entry_id=result.submission.inchikey,
                            name=name,
                            database_name=DATABASE_NAMES[source],
                            url=url,
                            raw=result.submission.smiles,
                            entry_type="compound",
                            primary_sequence=names,
                            fingerprint=fp,
                        )
                        if source == "mibig":
                            _apply_mibig_annotations(db, result.submission.inchikey, props, annotations, taxdb)
                        else:
                            type_label, genus, species = _npatlas_phylogeny(props)
                            resolution = resolve_phylogeny(
                                taxdb, genus=genus, species=species, fallback_type_label=type_label
                            )
                            db.add_phylogeny_annotation(
                                result.submission.inchikey,
                                type_label=resolution.type_label,
                                type_taxid=resolution.type_taxid,
                                genus=resolution.genus,
                                genus_taxid=resolution.genus_taxid,
                                species=resolution.species,
                                species_taxid=resolution.species_taxid,
                            )
                        added += 1

                    compounds += 1
                    pbar.update(1)
                    pbar.set_postfix(added=added, skipped=skipped)

                    if log_every > 0 and compounds % log_every == 0:
                        log.info(
                            "load_compounds[%s]: processed %d compounds (added=%d skipped=%d)",
                            source, compounds, added, skipped,
                        )
    finally:
        db.close()

    log.info("load_compounds[%s]: compounds=%d added=%d skipped=%d", source, compounds, added, skipped)


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
    ap.add_argument("--mibig-annotations", default=None, help="required when --source=mibig")
    ap.add_argument("--taxdump-dir", default=None, help="dir with NCBI names.dmp/nodes.dmp (skips taxid resolution if omitted)")
    ap.add_argument("--log-every", type=int, default=1000, help="log a progress line every N compounds (0 to disable)")
    args = ap.parse_args()

    run(
        results_path=args.results,
        db_path=args.db_path,
        source=args.source,
        reaction_rules_path=args.rxn_rules,
        matching_rules_path=args.mxn_rules,
        match_stereochemistry=args.match_stereochemistry,
        mibig_versions_path=args.mibig_versions,
        mibig_annotations_path=args.mibig_annotations,
        taxdump_dir=args.taxdump_dir,
        log_every=args.log_every,
    )


if __name__ == "__main__":
    main()
