"""Step 7: parse antiSMASH GenBank files into linear module readouts.

For each region: PARAS predicts NRPS A-domain substrate specificities
(retromol_antismash.inference.registry.annotate_region), then
retromol_antismash.modules.linear_readout collects PKS/NRPS modules in
biosynthetic order. One file per worker process (PARAS model loading is the
expensive per-process setup cost, so it's paid once per worker, not once per file).

PARAS is optional, gated by `paras_training_data_path`:
- None (default): PARAS is skipped entirely -- no model is loaded, annotate_region
  is never called -- and every NRPS module's substrate is set to an explicit
  "unknown" (see _mark_nrps_unknown) instead of a predicted one.
- set: raises NotImplementedError. On-the-fly retraining from a training data file
  isn't built yet; this exists so config can already carry the path for later.

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
import itertools
import json
import logging
from concurrent.futures import FIRST_COMPLETED, ProcessPoolExecutor, wait
from pathlib import Path

from rdkit import RDLogger
from tqdm import tqdm

from common import split_accession_version
from retromol_antismash.inference.model_paras import ParasModel
from retromol_antismash.inference.registry import annotate_region
from retromol_antismash.io import AntiSmashOptions, parse_antismash_gbk
from retromol_antismash.modules import LinearReadout, ModuleType, NRPSSubstrate, linear_readout

log = logging.getLogger(__name__)

GBK_GLOBS = ("*.gbk", "*.gb", "*.gbff")

UNKNOWN_NRPS_SUBSTRATE = NRPSSubstrate(name="unknown", smiles=None, score=0.0)

_G_PARAS_MODEL: ParasModel | None = None


def _mark_nrps_unknown(readout: LinearReadout) -> None:
    """Set every NRPS module's substrate to an explicit "unknown" -- used when PARAS is skipped
    so a module with no substrate call is visibly "not predicted", not indistinguishable from one
    PARAS looked at and genuinely couldn't call (module.predicted_substrate is mutable, not frozen)."""
    for module in readout.modules:
        if module.type == ModuleType.NRPS:
            module.predicted_substrate = UNKNOWN_NRPS_SUBSTRATE


def _init_worker(
    use_paras: bool, paras_threshold: float, paras_keep_top: int, paras_model_path: str | None, paras_cache_dir: str
) -> None:
    global _G_PARAS_MODEL
    # Belt-and-braces: importing this module (to resolve _init_worker as this pool's
    # initializer) already re-runs common.py's own RDLogger.DisableLog at the top of
    # this file's import chain, but this initializer is what the pool guarantees runs
    # once per worker no matter what -- see common.run_retromol_stream_quiet's
    # docstring for why that guarantee matters more than it might seem.
    RDLogger.DisableLog("rdApp.*")
    if use_paras:
        _G_PARAS_MODEL = ParasModel(
            threshold=paras_threshold,
            keep_top=paras_keep_top,
            model_path=paras_model_path,
            cache_dir=paras_cache_dir,
        )


def _process_file(path_str: str) -> tuple[list[dict], str | None]:
    """Parse one GenBank file. Returns (entries, error message)."""
    path = Path(path_str)
    try:
        raw_gbk = path.read_text()
        regions = parse_antismash_gbk(path, AntiSmashOptions())

        entries: list[dict] = []

        for region in regions:
            if _G_PARAS_MODEL is not None:
                annotate_region(region, domain_models=[_G_PARAS_MODEL])

            readout = linear_readout(region)

            if _G_PARAS_MODEL is None:
                _mark_nrps_unknown(readout)

            accession, _version = split_accession_version(region.id)

            entries.append({
                "accession": accession,
                "file_name": region.file_name,
                "raw_gbk": raw_gbk,
                "readout": readout.to_dict(),
            })

        return entries, None
    except Exception as e:
        return [], f"{path}: {e}"


def run(
    gbk_dir: str | Path,
    readouts_output_path: str | Path,
    paras_threshold: float = 0.1,
    paras_keep_top: int = 3,
    paras_model_path: str | Path | None = None,
    paras_cache_dir: str | Path = "paras_cache",
    paras_training_data_path: str | Path | None = None,
    workers: int = 1,
) -> None:
    if paras_training_data_path is not None:
        raise NotImplementedError(
            "retraining PARAS on the fly from a training data file is not implemented yet "
            f"(got paras_training_data_path={paras_training_data_path!r}); "
            "set paras.training_data_path back to null in config.yaml to skip PARAS "
            "and label every NRPS module 'unknown' instead"
        )

    gbk_dir = Path(gbk_dir)
    readouts_output_path = Path(readouts_output_path)
    readouts_output_path.parent.mkdir(parents=True, exist_ok=True)

    paths = sorted({p for pattern in GBK_GLOBS for p in gbk_dir.rglob(pattern)})

    n_files = 0
    n_entries = 0
    n_errors = 0

    # PARAS is always disabled for now -- paras_training_data_path is the only way
    # to opt in once retraining exists, and that path always raises above.
    use_paras = False

    init_args = (
        use_paras,
        paras_threshold,
        paras_keep_top,
        str(paras_model_path) if paras_model_path else None,
        str(paras_cache_dir),
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
    ap.add_argument("--paras-threshold", type=float, default=0.1)
    ap.add_argument("--paras-keep-top", type=int, default=3)
    ap.add_argument("--paras-model-path", default=None)
    ap.add_argument("--paras-cache-dir", default="paras_cache")
    ap.add_argument(
        "--paras-training-data-path",
        default=None,
        help="retrain PARAS on the fly from this file before parsing (NOT YET IMPLEMENTED -- raises if set)",
    )
    ap.add_argument("--workers", type=int, default=1)
    args = ap.parse_args()

    run(
        gbk_dir=args.gbk_dir,
        readouts_output_path=args.readouts_output,
        paras_threshold=args.paras_threshold,
        paras_keep_top=args.paras_keep_top,
        paras_model_path=args.paras_model_path,
        paras_cache_dir=args.paras_cache_dir,
        paras_training_data_path=args.paras_training_data_path,
        workers=args.workers,
    )


if __name__ == "__main__":
    main()
