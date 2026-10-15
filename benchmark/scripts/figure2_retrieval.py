"""Figure 2: top-1 retrieved compound's Tc for every MIBiG BGC.

For every MIBiG BGC (one query per accession; a BGC with multiple product compounds
is scored against whichever one the retrieval method itself ranks best, see
config's figure2.method), retrieve the single top-ranked compound candidate and
compute its whole-molecule Tanimoto (Tc) against the best-matching associated
ground-truth compound. Plotted as one density
per group: all MIBiG BGCs, plus each biosynthetic-class group configured in
benchmark/config.yaml's mibig_highlight_classes (default: Polyketide, Nonribosomal
peptide).
"""

import argparse
import json
import logging
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from tqdm import tqdm

from common import build_discovery_context, figure_style, normalize_for_alignment, open_db
from retrieval_engine import (
    TanimotoCache,
    build_eval_bgcs,
    combined_rank,
    fingerprint_only_rank,
    group_accessions_from_config,
    load_compound_universe,
    sample_candidate_pool,
    self_alignment_score,
    alignment_only_rank,
)

log = logging.getLogger(__name__)
_COLORS = ["#808080", "#4C72B0", "#DD8452", "#55A868"]


def _top1_tc(method: str, bq, ctx, db, universe, tanimoto: TanimotoCache, retrieval_cfg: dict, entries_by_id: dict) -> float | None:
    query_seq_norm = [normalize_for_alignment(n, ctx) for n in bq.primary_sequence]

    if method == "fingerprint_only":
        ranked = fingerprint_only_rank(db, bq.fingerprint, limit=1)
    elif method == "alignment_only":
        self_score = self_alignment_score(query_seq_norm, ctx)
        pool = sample_candidate_pool(
            universe,
            retrieval_cfg["alignment_only_candidate_pool_size"],
            retrieval_cfg["alignment_only_candidate_pool_seed"],
            must_include_ids=set(bq.ground_truth_ids),
        )
        ranked = alignment_only_rank(ctx, query_seq_norm, self_score, pool, retrieval_cfg["score_mode"], top_k=1)
    else:  # combined
        self_score = self_alignment_score(query_seq_norm, ctx)
        ranked = combined_rank(
            db, ctx, bq.fingerprint, query_seq_norm, self_score, retrieval_cfg["score_mode"],
            retrieval_cfg["fingerprint_prefilter_n"], top_k=1,
        )

    if not ranked:
        return None

    top1_smiles = entries_by_id.get(ranked[0].entry_id)
    if top1_smiles is None:
        return None

    best_tc = None
    for gt_id in bq.ground_truth_ids:
        gt_smiles = entries_by_id.get(gt_id)
        tc = tanimoto.tc(top1_smiles, gt_smiles)
        if tc is not None and (best_tc is None or tc > best_tc):
            best_tc = tc
    return best_tc


def run(
    db_path: str,
    links_path: str,
    method: str,
    mibig_highlight_classes: dict,
    retrieval_cfg: dict,
    tanimoto_cfg: dict,
    output_path: str,
) -> None:
    figure_style()

    with open(links_path) as fh:
        links = json.load(fh)

    ctx = build_discovery_context()
    universe = load_compound_universe(db_path)
    entries_by_id = {c.entry_id: c.smiles for c in universe}

    db = open_db(db_path)
    try:
        groups = build_eval_bgcs(db, links, group_accessions_from_config(db, links, mibig_highlight_classes))

        tanimoto = TanimotoCache(tanimoto_cfg["radius"], tanimoto_cfg["num_bits"], tanimoto_cfg["use_chirality"])

        tc_by_group: dict[str, list[float]] = {}
        for group_name, queries in groups.items():
            values = []
            for bq in tqdm(queries, desc=f"figure2[{group_name}]", unit="bgc"):
                tc = _top1_tc(method, bq, ctx, db, universe, tanimoto, retrieval_cfg, entries_by_id)
                if tc is not None:
                    values.append(tc)
            tc_by_group[group_name] = values
    finally:
        db.close()

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    x = np.linspace(0.0, 1.0, 400)
    labels = {"all": "All MIBiG BGCs", **{k: v["label"] for k, v in mibig_highlight_classes.items()}}
    for i, group_name in enumerate(["all", *mibig_highlight_classes.keys()]):
        values = tc_by_group.get(group_name, [])
        label = f"{labels[group_name]} (n={len(values)})"
        color = _COLORS[i % len(_COLORS)]
        if len(values) < 2 or np.std(values) == 0:
            if values:
                ax.axvline(values[0], color=color, label=label, linestyle="--" if group_name == "all" else "-")
            continue
        density = gaussian_kde(values)(x)
        ax.plot(x, density, color=color, linestyle="--" if group_name == "all" else "-", label=label)
        ax.fill_between(x, density, color=color, alpha=0.12)

    ax.set_xlim(0, 1)
    ax.set_xlabel("Top-1 retrieved compound Tc (vs. best-matching ground truth)")
    ax.set_ylabel("Density")
    ax.set_title(f"Figure 2 - Top-1 retrieval Tc per MIBiG BGC ({method})")
    ax.legend(loc="upper left", fontsize=8)
    fig.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)

    with open(output_path.with_suffix(".json"), "w") as fh:
        json.dump(tc_by_group, fh)

    log.info("figure2_retrieval: wrote %s", output_path)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--links", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    import yaml

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    run(
        db_path=args.db_path,
        links_path=args.links,
        method=cfg["figure2"]["method"],
        mibig_highlight_classes=cfg["mibig_highlight_classes"],
        retrieval_cfg=cfg["retrieval"],
        tanimoto_cfg=cfg["tanimoto"],
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
