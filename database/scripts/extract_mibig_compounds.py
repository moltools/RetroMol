"""Step 5 (prep): flatten MIBiG's per-BGC JSON files into one compounds JSONL.

MIBiG's JSON schema has shifted over releases (older releases nest everything
under a "cluster" key; newer ones are flatter), so field lookup here is
deliberately tolerant -- entries that don't have a resolvable accession/SMILES
are skipped (and counted) rather than raising.
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Any, Iterator

log = logging.getLogger(__name__)


def _entry_root(data: dict[str, Any]) -> dict[str, Any]:
    """Older MIBiG JSON nests everything under "cluster"; newer releases are flat."""
    cluster = data.get("cluster")
    return cluster if isinstance(cluster, dict) else data


def _accession(data: dict[str, Any], root: dict[str, Any]) -> str | None:
    return root.get("mibig_accession") or root.get("accession") or data.get("accession")


def _iter_compound_records(path: Path) -> Iterator[dict[str, Any]]:
    with open(path) as fh:
        data = json.load(fh)

    if not isinstance(data, dict):
        return

    root = _entry_root(data)
    accession = _accession(data, root)
    compounds = root.get("compounds")

    if not accession or not isinstance(compounds, list):
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


def run(mibig_json_dir: str | Path, output_path: str | Path) -> None:
    mibig_json_dir = Path(mibig_json_dir)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    json_files = sorted(mibig_json_dir.rglob("*.json"))
    written = 0
    skipped = 0

    with open(output_path, "w") as out:
        for path in json_files:
            try:
                had_any = False
                for record in _iter_compound_records(path):
                    out.write(json.dumps(record) + "\n")
                    written += 1
                    had_any = True
                if not had_any:
                    skipped += 1
            except Exception:
                log.exception("failed to parse MIBiG JSON file: %s", path)
                skipped += 1

    log.info("extract_mibig_compounds: wrote %d compound records, skipped %d files", written, skipped)


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mibig-json-dir", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    run(mibig_json_dir=args.mibig_json_dir, output_path=args.output)


if __name__ == "__main__":
    main()
