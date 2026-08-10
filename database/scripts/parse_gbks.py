"""Step 7: parse antiSMASH GenBank files into linear module readouts.

For each region: PARAS predicts NRPS A-domain substrate specificities
(retromol_antismash.inference.registry.annotate_region), then
retromol_antismash.modules.linear_readout collects PKS/NRPS modules in
biosynthetic order. One file per worker process (PARAS model loading is the
expensive per-process setup cost, so it's paid once per worker, not once per file).

Emits two outputs:
- readouts JSONL: one line per antiSMASH region, with its LinearReadout, the raw
  GenBank text of its source file, and the MIBiG accession/version/URL parsed out
  of the region id (e.g. "BGC0000001.5").
- a small accession -> version JSON map, reused by load_compounds.py to link MIBiG
  compound entries to the same BGC page.
"""

import argparse
import json
import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from common import mibig_url, split_accession_version
from retromol_antismash.inference.model_paras import ParasModel
from retromol_antismash.inference.registry import annotate_region
from retromol_antismash.io import AntiSmashOptions, parse_antismash_gbk
from retromol_antismash.modules import linear_readout

log = logging.getLogger(__name__)

GBK_GLOBS = ("*.gbk", "*.gb", "*.gbff")

_G_PARAS_MODEL: ParasModel | None = None


def _init_worker(paras_threshold: float, paras_keep_top: int, paras_model_path: str | None, paras_cache_dir: str) -> None:
    global _G_PARAS_MODEL
    _G_PARAS_MODEL = ParasModel(
        threshold=paras_threshold,
        keep_top=paras_keep_top,
        model_path=paras_model_path,
        cache_dir=paras_cache_dir,
    )


def _process_file(path_str: str) -> tuple[list[dict], dict[str, str], str | None]:
    """Parse one GenBank file. Returns (entries, accession->version, error message)."""
    path = Path(path_str)
    try:
        raw_gbk = path.read_text()
        regions = parse_antismash_gbk(path, AntiSmashOptions())

        entries: list[dict] = []
        versions: dict[str, str] = {}

        for region in regions:
            annotate_region(region, domain_models=[_G_PARAS_MODEL])
            readout = linear_readout(region)

            accession, version = split_accession_version(region.id)
            if version is not None:
                versions[accession] = version

            entries.append({
                "accession": accession,
                "version": version,
                "url": mibig_url(accession, version),
                "file_name": region.file_name,
                "raw_gbk": raw_gbk,
                "readout": readout.to_dict(),
            })

        return entries, versions, None
    except Exception as e:
        return [], {}, f"{path}: {e}"


def run(
    gbk_dir: str | Path,
    readouts_output_path: str | Path,
    versions_output_path: str | Path,
    paras_threshold: float = 0.1,
    paras_keep_top: int = 3,
    paras_model_path: str | Path | None = None,
    paras_cache_dir: str | Path = "paras_cache",
    workers: int = 1,
) -> None:
    gbk_dir = Path(gbk_dir)
    readouts_output_path = Path(readouts_output_path)
    versions_output_path = Path(versions_output_path)
    readouts_output_path.parent.mkdir(parents=True, exist_ok=True)
    versions_output_path.parent.mkdir(parents=True, exist_ok=True)

    paths = sorted({p for pattern in GBK_GLOBS for p in gbk_dir.rglob(pattern)})

    all_versions: dict[str, str] = {}
    n_entries = 0
    n_errors = 0

    init_args = (
        paras_threshold,
        paras_keep_top,
        str(paras_model_path) if paras_model_path else None,
        str(paras_cache_dir),
    )

    with open(readouts_output_path, "w", buffering=1) as out:
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=init_args) as pool:
            futures = [pool.submit(_process_file, str(p)) for p in paths]
            for fut in as_completed(futures):
                entries, versions, error = fut.result()
                if error is not None:
                    log.error("parse_gbks: %s", error)
                    n_errors += 1
                    continue
                for entry in entries:
                    out.write(json.dumps(entry) + "\n")
                    n_entries += 1
                all_versions.update(versions)

    with open(versions_output_path, "w") as fh:
        json.dump(all_versions, fh, indent=2, sort_keys=True)

    log.info("parse_gbks: files=%d regions=%d errors=%d", len(paths), n_entries, n_errors)


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gbk-dir", required=True)
    ap.add_argument("--readouts-output", required=True)
    ap.add_argument("--versions-output", required=True)
    ap.add_argument("--paras-threshold", type=float, default=0.1)
    ap.add_argument("--paras-keep-top", type=int, default=3)
    ap.add_argument("--paras-model-path", default=None)
    ap.add_argument("--paras-cache-dir", default="paras_cache")
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    run(
        gbk_dir=args.gbk_dir,
        readouts_output_path=args.readouts_output,
        versions_output_path=args.versions_output,
        paras_threshold=args.paras_threshold,
        paras_keep_top=args.paras_keep_top,
        paras_model_path=args.paras_model_path,
        paras_cache_dir=args.paras_cache_dir,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
