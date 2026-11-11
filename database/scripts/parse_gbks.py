"""Step 7: parse antiSMASH GenBank files into linear module readouts.

For each region: whichever model pmp.yml's "predictors.nrps.a_domain" currently
selects predicts NRPS A-domain substrate specificities
(retromol_antismash.inference.factory.build_nrps_a_domain_model +
retromol_antismash.inference.registry.annotate_region), then
retromol_antismash.modules.linear_readout collects PKS/NRPS modules in
biosynthetic order -- using the SAME factory + pmp.yml the GUI backend
(gui/src/server/routes/jobs.py) and the CLI test script (scripts/predict_bgc.py)
use, so all three ways of parsing a GenBank file resolve substrates identically.
One file per worker process (model loading is the expensive per-process setup
cost, so it's paid once per worker, not once per file).

If pmp.yml's selected method is `source: qualifier` instead of `source: model`
(e.g. reading antiSMASH's own NRPS substrate call straight from the GenBank
record), no model is built or run at all -- collect_nrps_modules reads the
qualifier directly, same as everywhere else this config is used.

Emits one output: readouts JSONL, one line per antiSMASH region, with its
LinearReadout, the raw GenBank text of its source file, and the MIBiG accession
parsed out of the region id.

MIBiG's URL (built from an accession + a per-entry revision number, e.g.
".../BGC0000001.5/...") used to be derivable right here, from the GBK's own
ACCESSION.VERSION line -- but MIBiG 4.0's GBKs dropped the version suffix
entirely ("ACCESSION   BGC0000001", no VERSION line). The version now only
exists in the JSON's own top-level "version" field, so extract_mibig_compounds.py
is the sole source of the accession -> version map (mibig_versions.json) both
this script's own consumer (load_bgcs.py) and load_compounds.py need -- this
script only emits `accession` (still correct: it's just region.id, with or
without a version suffix) for load_bgcs.py to look up.
"""

import argparse
import hashlib
import itertools
import json
import logging
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

from rdkit import RDLogger
from tqdm import tqdm

from common import split_accession_version
from retromol_antismash.inference.base import DomainInferenceModel
from retromol_antismash.inference.factory import build_nrps_a_domain_model
from retromol_antismash.inference.registry import annotate_region
from retromol_antismash.io import AntiSmashOptions, parse_antismash_gbk
from retromol_antismash.modules import linear_readout
from retromol_antismash.predictions import PredictionConfig

log = logging.getLogger(__name__)

GBK_GLOBS = ("*.gbk", "*.gb", "*.gbff")

_G_CONFIG: PredictionConfig | None = None
_G_NRPS_MODEL: DomainInferenceModel | None = None


def _init_worker(
    pmp_path: str | None,
    threshold: float,
    keep_top: int,
    cache_dir: str,
    force_retrain: bool,
) -> None:
    global _G_CONFIG, _G_NRPS_MODEL
    # Belt-and-braces: importing this module (to resolve _init_worker as this pool's
    # initializer) already re-runs common.py's own RDLogger.DisableLog at the top of
    # this file's import chain, but this initializer is what the pool guarantees runs
    # once per worker no matter what -- see common.run_retromol_stream_quiet's
    # docstring for why that guarantee matters more than it might seem.
    RDLogger.DisableLog("rdApp.*")
    _G_CONFIG = PredictionConfig.load_from_file(pmp_path) if pmp_path else PredictionConfig.load_default()
    # None when pmp.yml's predictors.nrps.a_domain is source: qualifier -- nothing to
    # build or register, collect_nrps_modules reads the qualifier directly.
    _G_NRPS_MODEL = build_nrps_a_domain_model(
        _G_CONFIG,
        threshold=threshold,
        keep_top=keep_top,
        cache_dir=cache_dir,
        force_retrain=force_retrain,
    )


def _process_file(path_str: str) -> tuple[list[dict], str | None]:
    """Parse one GenBank file. Returns (entries, error message)."""
    path = Path(path_str)
    try:
        raw_gbk = path.read_text()
        file_hash = hashlib.sha256(raw_gbk.encode("utf-8")).hexdigest()
        regions = parse_antismash_gbk(path, AntiSmashOptions())

        entries: list[dict] = []

        for region in regions:
            if _G_NRPS_MODEL is not None:
                annotate_region(region, domain_models=[_G_NRPS_MODEL])

            readout = linear_readout(region, config=_G_CONFIG)

            accession, _version = split_accession_version(region.id)

            entries.append({
                "accession": accession,
                "file_name": region.file_name,
                "raw_gbk": raw_gbk,
                "file_hash": file_hash,
                "readout": readout.to_dict(),
            })

        return entries, None
    except Exception as e:
        return [], f"{path}: {e}"


def run(
    gbk_dir: str | Path,
    readouts_output_path: str | Path,
    pmp_path: str | Path | None = None,
    paras_threshold: float = 0.1,
    paras_keep_top: int = 3,
    paras_cache_dir: str | Path = "paras_cache",
    force_retrain: bool = False,
    workers: int = 1,
) -> None:
    """
    :param gbk_dir: directory of antiSMASH GenBank files to parse (searched recursively).
    :param readouts_output_path: where to write the readouts JSONL.
    :param pmp_path: path to a pmp.yml prediction-mapping file, or None to use the packaged
        default -- see retromol_antismash.inference.factory.build_nrps_a_domain_model for
        how this selects (or skips) an NRPS substrate-prediction model.
    :param paras_threshold: minimum predicted probability for a substrate call to be kept.
    :param paras_keep_top: number of top-scoring substrate predictions to keep per domain.
    :param paras_cache_dir: training-signature + fitted-model cache directory (see
        retromol_paras.train.train_model).
    :param force_retrain: retrain from scratch even if a cached model exists.
    :param workers: number of worker processes.
    """
    gbk_dir = Path(gbk_dir)
    readouts_output_path = Path(readouts_output_path)
    readouts_output_path.parent.mkdir(parents=True, exist_ok=True)

    paths = sorted({p for pattern in GBK_GLOBS for p in gbk_dir.rglob(pattern)})

    n_files = 0
    n_entries = 0
    n_errors = 0

    init_args = (
        str(pmp_path) if pmp_path else None,
        paras_threshold,
        paras_keep_top,
        str(paras_cache_dir),
        force_retrain,
    )

    # Bounded sliding window rather than submitting every file as a future up front --
    # at hundreds of thousands of files, submitting everything at once means that many
    # Future/WorkItem objects (plus each one's eventual result: raw GBK text + a full
    # serialized LinearReadout) sit in memory simultaneously. Keeping at most
    # `max_pending` in flight bounds memory to a small multiple of `workers`,
    # regardless of corpus size -- same shape as retromol.io.streaming's own batching
    # that parse_compounds.py relies on for the same reason.
    max_pending = max(workers * 4, 1)
    paths_iter = iter(paths)

    with open(readouts_output_path, "w", buffering=1) as out:
        with ProcessPoolExecutor(max_workers=workers, initializer=_init_worker, initargs=init_args) as pool:
            with tqdm(total=len(paths), desc="parse_gbks", unit="file") as pbar:
                pending = {pool.submit(_process_file, str(p)) for p in itertools.islice(paths_iter, max_pending)}

                while pending:
                    done, pending = wait(pending, return_when=FIRST_COMPLETED)

                    for fut in done:
                        entries, error = fut.result()
                        n_files += 1

                        if error is not None:
                            log.error("parse_gbks: %s", error)
                            n_errors += 1
                        else:
                            for entry in entries:
                                out.write(json.dumps(entry) + "\n")
                                n_entries += 1

                        pbar.update(1)
                        pbar.set_postfix(regions=n_entries, errors=n_errors)

                        next_path = next(paths_iter, None)
                        if next_path is not None:
                            pending.add(pool.submit(_process_file, str(next_path)))

    log.info("parse_gbks: files=%d regions=%d errors=%d", n_files, n_entries, n_errors)


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gbk-dir", required=True)
    ap.add_argument("--readouts-output", required=True)
    ap.add_argument("--pmp", default=None, help="path to a pmp.yml prediction-mapping file (default: packaged pmp.yml)")
    ap.add_argument("--paras-threshold", type=float, default=0.1)
    ap.add_argument("--paras-keep-top", type=int, default=3)
    ap.add_argument("--paras-cache-dir", default="paras_cache")
    ap.add_argument("--force-retrain", action="store_true", help="retrain from scratch even if a cached model exists")
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    run(
        gbk_dir=args.gbk_dir,
        readouts_output_path=args.readouts_output,
        pmp_path=args.pmp,
        paras_threshold=args.paras_threshold,
        paras_keep_top=args.paras_keep_top,
        paras_cache_dir=args.paras_cache_dir,
        force_retrain=args.force_retrain,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
