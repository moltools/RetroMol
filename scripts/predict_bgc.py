#!/usr/bin/env python3
"""
Parse an antiSMASH-annotated GenBank file and print the PKS/NRPS substrate
predictions RetroMol resolves for each module -- a quick way to check what
`pmp.yml` (see src/retromol/data/pmp.yml) actually produces against a real
GBK without going through the full compound-matching/alignment pipeline.

Examples:

    # Basic run, using the packaged default pmp.yml (PKS from antiSMASH
    # qualifiers, NRPS from PARAS):
    python scripts/predict_bgc.py path/to/cluster.region001.gbk

    # Skip the NRPS model entirely (e.g. no NRPS in the cluster, or you just want
    # antiSMASH's own qualifier fallback without downloading/running/training a model):
    python scripts/predict_bgc.py path/to/cluster.gbk --no-nrps-model

    # Point at an already-trained model cache instead of training fresh:
    python scripts/predict_bgc.py path/to/cluster.gbk \\
        --paras-cache-dir /path/to/paras_cache

    # Try edits to pmp.yml/mxn.yml before committing them:
    python scripts/predict_bgc.py path/to/cluster.gbk \\
        --pmp path/to/my_pmp.yml --mxn path/to/my_mxn.yml

    # Machine-readable output:
    python scripts/predict_bgc.py path/to/cluster.gbk --json out.json
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from pathlib import Path

from retromol.model.rules import RuleSet
from retromol_antismash.inference.factory import build_nrps_a_domain_model
from retromol_antismash.inference.registry import annotate_region, register_domain_model
from retromol_antismash.io import AntiSmashOptions, parse_antismash_gbk
from retromol_antismash.model import Region
from retromol_antismash.modules import (
    LinearReadout,
    Module,
    ModuleType,
    linear_readout,
    module_primary_sequence_tokens,
)
from retromol_antismash.predictions import PredictionConfig


log = logging.getLogger("predict_bgc")


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Print RetroMol's PKS/NRPS substrate predictions for an antiSMASH GenBank file.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument("gbk", type=Path, help="path to an antiSMASH-annotated GenBank file")
    parser.add_argument(
        "--readout-level",
        choices=["region", "cand_cluster"],
        default="region",
        help="antiSMASH feature type to treat as a biosynthetic cluster (default: region)",
    )

    rules = parser.add_argument_group("rule files")
    rules.add_argument("--pmp", type=Path, default=None, help="path to a pmp.yml prediction-mapping file (default: packaged pmp.yml)")
    rules.add_argument("--mxn", type=Path, default=None, help="path to a mxn.yml matching-rules file (default: packaged mxn.yml)")
    rules.add_argument("--rxn", type=Path, default=None, help="path to a rxn.yml reaction-rules file (default: packaged rxn.yml)")

    nrps = parser.add_argument_group("NRPS model (used when pmp.yml selects a source: model method for nrps.a_domain, e.g. paras/paras_cli)")
    nrps.add_argument("--no-nrps-model", "--no-paras", dest="no_nrps_model", action="store_true", help="don't register/run any NRPS model, even if pmp.yml selects one")
    nrps.add_argument("--paras-cache-dir", type=Path, default=None, help="training-signature + fitted-model cache directory for 'paras_cli' (see retromol_paras.train.train_model)")
    nrps.add_argument("--paras-threshold", type=float, default=0.1, help="probability threshold (default: 0.1)")
    nrps.add_argument("--paras-keep-top", type=int, default=3, help="number of top predictions to keep per domain (default: 3)")
    nrps.add_argument("--force-retrain", action="store_true", help="'paras_cli'-only: retrain from scratch even if a cached model exists")

    out = parser.add_argument_group("output")
    out.add_argument("--json", type=Path, default=None, metavar="PATH", help="also write machine-readable results to this JSON file")
    out.add_argument("-l", "--log-level", choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"], default="WARNING")

    return parser


def _maybe_register_nrps_model(config: PredictionConfig, args: argparse.Namespace) -> None:
    """
    Register whichever NRPS A-domain model pmp.yml's "predictors.nrps.a_domain"
    currently selects (unless --no-nrps-model was passed) -- the same factory
    the database pipeline (database/scripts/parse_gbks.py) and the GUI backend
    (gui/src/server/routes/jobs.py) use, so this script resolves substrates
    identically to both.
    """
    if args.no_nrps_model:
        return

    model = build_nrps_a_domain_model(
        config,
        threshold=args.paras_threshold,
        keep_top=args.paras_keep_top,
        cache_dir=args.paras_cache_dir,
        force_retrain=args.force_retrain,
    )
    if model is None:
        return

    log.info(f"registering {model.name!r} model for NRPS A-domain substrate prediction")
    register_domain_model(model)


def _module_report(module: Module, ruleset: RuleSet) -> dict:
    name, tokens = module_primary_sequence_tokens(module, ruleset)
    row: dict = {
        "gene_id": module.gene_id,
        "module_index_in_gene": module.module_index_in_gene,
        "start": module.start,
        "end": module.end,
        "type": module.type.value,
        "present_domains": module.present_domains,
        "resolved_name": name,
        "resolved_tokens": tokens,
    }

    if module.type is ModuleType.PKS:
        substrate = module.substrate
        row["anatomy"] = module.anatomy.to_dict()
        row["extender_unit"] = substrate.extender_unit.value
        row["substituent_digit"] = substrate.substituent_type
    else:
        substrate = module.substrate
        row["anatomy"] = module.anatomy.to_dict()
        row["predicted_substrate"] = substrate.to_dict() if substrate else None

    return row


def _print_report(readout: LinearReadout, rows: list[dict]) -> None:
    print(f"\n=== {readout.id} ({readout.file_name}) [{readout.start}-{readout.end}] ===")
    if readout.modifiers:
        print(f"modifiers: {', '.join(readout.modifiers)}")

    if not rows:
        print("(no PKS/NRPS modules found)")
        return

    for row in rows:
        loc = f"{row['gene_id']}[{row['module_index_in_gene']}]"
        if row["type"] == ModuleType.PKS.value:
            letter = row["extender_unit"].rsplit("_", 1)[-1]  # "PK_B" -> "B"
            digit = row["substituent_digit"] or "?"
            detail = f"extender={letter}{digit} ({row['extender_unit']})"
        else:
            sub = row["predicted_substrate"]
            if sub:
                detail = f"substrate={sub['name']!r} score={sub['score']}"
            else:
                detail = "substrate=(none)"
        print(f"  {loc:<28} {row['type']:<4} {detail:<40} -> {row['resolved_name']:<10} tokens={row['resolved_tokens']}")


def main() -> None:
    args = build_arg_parser().parse_args()
    logging.basicConfig(level=args.log_level, format="%(levelname)s %(name)s: %(message)s")

    if not args.gbk.exists():
        log.error(f"no such file: {args.gbk}")
        sys.exit(1)

    config = PredictionConfig.load_from_file(args.pmp) if args.pmp else PredictionConfig.load_default()
    ruleset = RuleSet.load_from_files(reaction_rules_path=args.rxn, matching_rules_path=args.mxn)

    _maybe_register_nrps_model(config, args)

    options = AntiSmashOptions(readout_level=args.readout_level)
    regions: list[Region] = parse_antismash_gbk(args.gbk, options)
    if not regions:
        log.warning(f"no {args.readout_level!r}-level features found in {args.gbk}")

    all_rows: list[dict] = []
    for region in regions:
        annotate_region(region)  # runs any registered gene/domain inference models (e.g. PARAS)
        readout = linear_readout(region, config=config)

        rows = [_module_report(m, ruleset) for m in readout.biosynthetic_order()]
        _print_report(readout, rows)

        all_rows.append({
            "region_id": readout.id,
            "file_name": readout.file_name,
            "start": readout.start,
            "end": readout.end,
            "modifiers": readout.modifiers,
            "modules": rows,
        })

    if args.json:
        args.json.write_text(json.dumps(all_rows, indent=2, default=str))
        log.info(f"wrote {args.json}")


if __name__ == "__main__":
    main()
