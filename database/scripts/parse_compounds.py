"""Steps 3 & 5: run RetroMol over a batch of compounds (NPAtlas SDF or MIBiG compounds JSONL).

Shared between both compound sources -- only the input format differs. This is the
slow step (one retrosynthetic analysis per compound), so it's parallelized with
`workers` worker processes via common.run_retromol_stream_quiet -- the same
worker/batching machinery the `retromol` CLI's batch mode uses
(retromol.io.streaming.run_retromol_stream), but with RDKit's C-level logging also
disabled inside every worker (see common.py for why that needs its own initializer).
"""

import argparse
import json
import logging
from pathlib import Path
from typing import Literal

from tqdm import tqdm

from common import load_ruleset, run_retromol_stream_quiet
from retromol.io.streaming import stream_json_records, stream_sdf_records

log = logging.getLogger(__name__)


def run(
    input_path: str | Path,
    input_format: Literal["sdf", "jsonl"],
    output_path: str | Path,
    reaction_rules_path: str | Path | None,
    matching_rules_path: str | Path | None,
    match_stereochemistry: bool = False,
    smiles_col: str = "smiles",
    workers: int = 1,
    batch_size: int = 2000,
    log_every: int = 1000,
) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    ruleset = load_ruleset(reaction_rules_path, matching_rules_path, match_stereochemistry)

    if input_format == "sdf":
        row_iter = stream_sdf_records(str(input_path))
    else:
        row_iter = stream_json_records(str(input_path), jsonl=True)

    successes = 0
    errors = 0

    with open(output_path, "w", buffering=1) as out:
        with tqdm(desc="parse_compounds", unit="cmpd") as pbar:
            for evt in run_retromol_stream_quiet(
                ruleset=ruleset,
                row_iter=row_iter,
                smiles_col=smiles_col,
                workers=workers,
                batch_size=batch_size,
            ):
                if evt.error is not None:
                    errors += 1
                elif evt.result is not None:
                    out.write(json.dumps(evt.result) + "\n")
                    successes += 1

                pbar.update(1)
                pbar.set_postfix(ok=successes, err=errors)

                total = successes + errors
                if log_every > 0 and total % log_every == 0:
                    log.info("parse_compounds: parsed %d (successes=%d errors=%d)", total, successes, errors)

    log.info("parse_compounds: successes=%d errors=%d -> %s", successes, errors, output_path)


def main() -> None:
    logging.basicConfig(level=logging.INFO)

    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True)
    ap.add_argument("--input-format", choices=["sdf", "jsonl"], required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--rxn-rules", default=None)
    ap.add_argument("--mxn-rules", default=None)
    ap.add_argument("--match-stereochemistry", action="store_true")
    ap.add_argument("--smiles-col", default="smiles")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--batch-size", type=int, default=2000)
    ap.add_argument("--log-every", type=int, default=1000, help="log a progress line every N compounds (0 to disable)")
    args = ap.parse_args()

    run(
        input_path=args.input,
        input_format=args.input_format,
        output_path=args.output,
        reaction_rules_path=args.rxn_rules,
        matching_rules_path=args.mxn_rules,
        match_stereochemistry=args.match_stereochemistry,
        smiles_col=args.smiles_col,
        workers=args.workers,
        batch_size=args.batch_size,
        log_every=args.log_every,
    )


if __name__ == "__main__":
    main()
