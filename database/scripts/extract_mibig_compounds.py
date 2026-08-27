"""Step 5 (prep): flatten MIBiG's per-BGC JSON files into one compounds JSONL,
plus an accession -> version map (also consumed by load_bgcs.py, see versions_output).

MIBiG's JSON schema has shifted over releases (older releases nest everything
under a "cluster" key; newer ones are flatter), so field lookup here is
deliberately tolerant -- entries that don't have a resolvable accession/SMILES
are skipped (and counted) rather than raising.

The accession -> version map used to come from the GBK files' own
ACCESSION.VERSION line (see common.split_accession_version) -- that broke with
MIBiG 4.0, whose GBKs dropped the version suffix entirely ("ACCESSION
BGC0000001", no "VERSION" line at all). The JSON files still carry it, as a
top-level "version" int (e.g. 5, matching the ".5" in
https://mibig.secondarymetabolites.org/repository/BGC0000001.5/index.html),
so that's the source of truth now -- extracted here, unconditionally (even for
files with no usable compound records), since parse_gbks.py can no longer
provide it.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Any, Iterator

from tqdm import tqdm

log = logging.getLogger(__name__)


def _entry_root(data: dict[str, Any]) -> dict[str, Any]:
    """Older MIBiG JSON nests everything under "cluster"; newer releases are flat."""
    cluster = data.get("cluster")
    return cluster if isinstance(cluster, dict) else data


def _accession(data: dict[str, Any], root: dict[str, Any]) -> str | None:
    return root.get("mibig_accession") or root.get("accession") or data.get("accession")


def _version(data: dict[str, Any], root: dict[str, Any]) -> str | None:
    version = data.get("version")
    if version is None:
        version = root.get("version")
    return str(version) if version is not None else None


def _iter_compound_records(data: dict[str, Any], root: dict[str, Any], accession: str) -> Iterator[dict[str, Any]]:
    compounds = root.get("compounds")
    if not isinstance(compounds, list):
        return

    for idx, compound in enumerate(compounds):
        if not isinstance(compound, dict):
            continue

        smiles = compound.get("chem_struct") or compound.get("smiles") or compound.get("structure")
        if not smiles:
            continue

        name = compound.get("compound") or compound.get("name") or f"{accession} compound {idx + 1}"

        yield {
            "id": f"{accession}:{idx}",
            "smiles": smiles,
            "name": name,
            "mibig_accession": accession,
        }


# MIBiG's own short biosynthesis-class codes -> friendlier display labels for the
# chemical_class annotation. Covers every value observed across the full 4.0 corpus
# (PKS/NRPS/ribosomal/other/terpene/saccharide); an unrecognized future code falls
# back to itself unchanged rather than being dropped.
BIOSYN_CLASS_LABELS = {
    "PKS": "Polyketide",
    "NRPS": "Nonribosomal peptide",
    "ribosomal": "RiPP",
    "terpene": "Terpene",
    "saccharide": "Saccharide",
    "other": "Other",
}


def _annotations(root: dict[str, Any]) -> dict[str, Any]:
    """Phylogeny + chemical-class metadata shared by every compound/BGC under one accession.

    MIBiG 4.0's real (flat) schema nests these under "taxonomy" ({"name", "ncbiTaxId"})
    and "biosynthesis" ({"classes": [{"class": "PKS", ...}, ...]}) -- not the top-level
    "organism_name"/"ncbi_tax_id"/"biosyn_class" keys older MIBiG releases used.
    """
    taxonomy = root.get("taxonomy")
    taxonomy = taxonomy if isinstance(taxonomy, dict) else {}
    organism_name = taxonomy.get("name")
    ncbi_tax_id = taxonomy.get("ncbiTaxId")

    biosynthesis = root.get("biosynthesis")
    biosynthesis = biosynthesis if isinstance(biosynthesis, dict) else {}
    classes = biosynthesis.get("classes")
    classes = classes if isinstance(classes, list) else []
    biosyn_class = [
        BIOSYN_CLASS_LABELS.get(c.get("class"), c.get("class"))
        for c in classes
        if isinstance(c, dict) and c.get("class")
    ]

    return {
        "organism_name": organism_name if isinstance(organism_name, str) else None,
        "ncbi_tax_id": str(ncbi_tax_id) if ncbi_tax_id is not None else None,
        "biosyn_class": biosyn_class,
    }


def run(
    mibig_json_dir: str | Path,
    output_path: str | Path,
    versions_output_path: str | Path,
    annotations_output_path: str | Path,
) -> None:
    mibig_json_dir = Path(mibig_json_dir)
    output_path = Path(output_path)
    versions_output_path = Path(versions_output_path)
    annotations_output_path = Path(annotations_output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    versions_output_path.parent.mkdir(parents=True, exist_ok=True)
    annotations_output_path.parent.mkdir(parents=True, exist_ok=True)

    json_files = sorted(mibig_json_dir.rglob("*.json"))
    written = 0
    skipped = 0
    versions: dict[str, str] = {}
    annotations: dict[str, dict[str, Any]] = {}

    with open(output_path, "w") as out:
        for path in tqdm(json_files, desc="extract_mibig_compounds", unit="file"):
            try:
                with open(path) as fh:
                    data = json.load(fh)

                if not isinstance(data, dict):
                    skipped += 1
                    continue

                root = _entry_root(data)
                accession = _accession(data, root)
                if not accession:
                    skipped += 1
                    continue

                version = _version(data, root)
                if version is not None:
                    versions[accession] = version

                annotations[accession] = _annotations(root)

                had_any = False
                for record in _iter_compound_records(data, root, accession):
                    out.write(json.dumps(record) + "\n")
                    written += 1
                    had_any = True
                if not had_any:
                    skipped += 1
            except Exception:
                log.exception("failed to parse MIBiG JSON file: %s", path)
                skipped += 1

    with open(versions_output_path, "w") as fh:
        json.dump(versions, fh, indent=2, sort_keys=True)

    with open(annotations_output_path, "w") as fh:
        json.dump(annotations, fh, indent=2, sort_keys=True)

    log.info(
        "extract_mibig_compounds: wrote %d compound records, skipped %d files, "
        "resolved %d accession versions, %d accession annotations",
        written, skipped, len(versions), len(annotations),
    )


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mibig-json-dir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--versions-output", required=True)
    ap.add_argument("--annotations-output", required=True)
    args = ap.parse_args()

    run(
        mibig_json_dir=args.mibig_json_dir,
        output_path=args.output,
        versions_output_path=args.versions_output,
        annotations_output_path=args.annotations_output,
    )


if __name__ == "__main__":
    main()
