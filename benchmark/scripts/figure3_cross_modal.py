"""Figure 3: cross-modal retrieval hit@k, and per-method retrieval timing.

For every MIBiG BGC query (grouped as in figure2: all, plus each biosynthetic-class
group in benchmark/config.yaml's mibig_highlight_classes), retrieve up to max(top_ks)
compound candidates with each of three methods (fingerprint-only, alignment-only,
combined fingerprint-prefilter-then-realign; see retrieval_engine.py), and check
whether the best candidate among the top-k has a whole-molecule Tc >= each configured
threshold against any of that BGC's ground-truth compounds. Hit rate is plotted vs.
k (one line per method) for every (group, threshold) pair. Wall-clock time per query
is also recorded per method and plotted separately (figure3_timing.png), measured
once, over the "all" group (retrieval cost per query doesn't depend on which
biosynthetic-class group it's reported under).
"""

import argparse
import json
import logging
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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
    timed_rank,
    alignment_only_rank,
)

log = logging.getLogger(__name__)
_METHOD_LABELS = {"fingerprint_only": "Fingerprint only", "alignment_only": "Pairwise alignment only", "combined": "Fingerprint + alignment"}
_METHOD_COLORS = {"fingerprint_only": "#4C72B0", "alignment_only": "#DD8452", "combined": "#55A868"}


def _rank_for_method(method: str, bq, ctx, db, universe, retrieval_cfg: dict, max_k: int):
    query_seq_norm = [normalize_for_alignment(n, ctx) for n in bq.primary_sequence]

    if method == "fingerprint_only":
        return timed_rank(fingerprint_only_rank, db, bq.fingerprint, max_k)
    if method == "alignment_only":
        self_score = self_alignment_score(query_seq_norm, ctx)
        pool = sample_candidate_pool(
            universe,
            retrieval_cfg["alignment_only_candidate_pool_size"],
            retrieval_cfg["alignment_only_candidate_pool_seed"],
            must_include_ids=set(bq.ground_truth_ids),
        )
        return timed_rank(alignment_only_rank, ctx, query_seq_norm, self_score, pool, retrieval_cfg["score_mode"], max_k)
    # combined
    self_score = self_alignment_score(query_seq_norm, ctx)
    return timed_rank(
        combined_rank, db, ctx, bq.fingerprint, query_seq_norm, self_score,
        retrieval_cfg["score_mode"], retrieval_cfg["fingerprint_prefilter_n"], max_k,
    )


def _plot_hit_rates(
    hit_rates: dict, group_labels: dict, methods: list, top_ks: list, thresholds: list, output_path: str
) -> None:
    groups = list(hit_rates.keys())
    fig, axes = plt.subplots(len(groups), len(thresholds), figsize=(5.5 * len(thresholds), 3.8 * len(groups)), squeeze=False)
    for row, group_name in enumerate(groups):
        for col, threshold in enumerate(thresholds):
            ax = axes[row][col]
            for method in methods:
                ys = [hit_rates[group_name][method][str(threshold)][str(k)] for k in top_ks]
                ax.plot(top_ks, ys, marker="o", color=_METHOD_COLORS[method], label=_METHOD_LABELS[method])
            ax.set_xticks(top_ks)
            ax.set_ylim(0, 1.05)
            ax.set_title(f"{group_labels[group_name]} - Tc >= {threshold}")
            ax.set_xlabel("k")
            if col == 0:
                ax.set_ylabel("Hit rate")
            if row == 0 and col == len(thresholds) - 1:
                ax.legend(loc="lower right", fontsize=8)
    fig.suptitle("Figure 3 - Cross-modal retrieval hit@k")
    fig.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def _plot_timing(timings: dict, methods: list, output_path: str) -> None:
    fig2, ax2 = plt.subplots(figsize=(5, 4))
    means = [float(np.mean(timings[m])) if timings[m] else 0.0 for m in methods]
    stds = [float(np.std(timings[m])) if timings[m] else 0.0 for m in methods]
    ax2.bar([_METHOD_LABELS[m] for m in methods], means, yerr=stds, color=[_METHOD_COLORS[m] for m in methods])
    ax2.set_ylabel("Mean retrieval time per BGC query (s)")
    ax2.set_title(f"Figure 3 - Retrieval timing (n={len(timings[methods[0]])} queries)")
    plt.setp(ax2.get_xticklabels(), rotation=20, ha="right")
    fig2.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig2.savefig(output_path)
    plt.close(fig2)


def run(
    db_path: str,
    links_path: str,
    figure3_cfg: dict,
    mibig_highlight_classes: dict,
    retrieval_cfg: dict,
    tanimoto_cfg: dict,
    output_path: str,
    timing_output_path: str,
) -> None:
    figure_style()

    with open(links_path) as fh:
        links = json.load(fh)

    group_labels = {"all": "All MIBiG BGCs", **{k: v["label"] for k, v in mibig_highlight_classes.items()}}

    ctx = build_discovery_context()
    universe = load_compound_universe(db_path)
    entries_by_id = {c.entry_id: c.smiles for c in universe}
    tanimoto = TanimotoCache(tanimoto_cfg["radius"], tanimoto_cfg["num_bits"], tanimoto_cfg["use_chirality"])

    top_ks = figure3_cfg["top_ks"]
    thresholds = figure3_cfg["tc_thresholds"]
    methods = figure3_cfg["methods"]
    max_k = max(top_ks)

    db = open_db(db_path)
    try:
        groups = build_eval_bgcs(db, links, group_accessions_from_config(db, links, mibig_highlight_classes))

        # hit_rates[group][method][threshold][k] = fraction of queries with a hit
        hit_rates: dict[str, dict[str, dict[float, dict[int, float]]]] = {}
        timings: dict[str, list[float]] = {m: [] for m in methods}

        for group_name, queries in groups.items():
            hit_rates[group_name] = {m: {t: {} for t in thresholds} for m in methods}
            for method in methods:
                per_query_hits: dict[int, dict[float, list[bool]]] = {k: {t: [] for t in thresholds} for k in top_ks}

                for bq in tqdm(queries, desc=f"figure3[{group_name}][{method}]", unit="bgc"):
                    ranked, elapsed = _rank_for_method(method, bq, ctx, db, universe, retrieval_cfg, max_k)
                    if group_name == "all":
                        timings[method].append(elapsed)

                    candidate_best_tc = []
                    for cand in ranked:
                        cand_smiles = entries_by_id.get(cand.entry_id)
                        best_tc = 0.0
                        for gt_id in bq.ground_truth_ids:
                            tc = tanimoto.tc(cand_smiles, entries_by_id.get(gt_id))
                            if tc is not None:
                                best_tc = max(best_tc, tc)
                        candidate_best_tc.append(best_tc)

                    for k in top_ks:
                        topk_best = max(candidate_best_tc[:k], default=0.0)
                        for threshold in thresholds:
                            per_query_hits[k][threshold].append(topk_best >= threshold)

                for k in top_ks:
                    for threshold in thresholds:
                        vals = per_query_hits[k][threshold]
                        hit_rates[group_name][method][threshold][k] = (
                            float(np.mean(vals)) if vals else float("nan")
                        )
    finally:
        db.close()

    # _plot_hit_rates always indexes by str(threshold)/str(k), matching how JSON
    # round-trips keys (see replot()) -- normalize here so both callers agree.
    hit_rates_json = {
        group_name: {
            method: {str(t): {str(k): v for k, v in kmap.items()} for t, kmap in tmap.items()}
            for method, tmap in mmap.items()
        }
        for group_name, mmap in hit_rates.items()
    }

    output_path = Path(output_path)
    timing_output_path = Path(timing_output_path)

    _plot_hit_rates(hit_rates_json, group_labels, methods, top_ks, thresholds, output_path)
    with open(output_path.with_suffix(".json"), "w") as fh:
        json.dump(hit_rates_json, fh, indent=2)

    _plot_timing(timings, methods, timing_output_path)
    with open(timing_output_path.with_suffix(".json"), "w") as fh:
        json.dump({m: timings[m] for m in methods}, fh)

    log.info("figure3_cross_modal: wrote %s and %s", output_path, timing_output_path)


def replot(
    hit_rates_json_path: str,
    timing_json_path: str,
    figure3_cfg: dict,
    mibig_highlight_classes: dict,
    output_path: str,
    timing_output_path: str,
) -> None:
    """Redraw both figures straight from a previous run's saved JSON (see run()'s
    output_path.with_suffix(".json") / timing_output_path.with_suffix(".json")) --
    for a plotting-only change (styling, axis scale, labels), so the expensive
    retrieval computation doesn't need to be redone just to redraw the same numbers.
    """
    figure_style()

    with open(hit_rates_json_path) as fh:
        hit_rates = json.load(fh)
    with open(timing_json_path) as fh:
        timings = json.load(fh)

    group_labels = {"all": "All MIBiG BGCs", **{k: v["label"] for k, v in mibig_highlight_classes.items()}}
    top_ks = figure3_cfg["top_ks"]
    thresholds = figure3_cfg["tc_thresholds"]
    methods = figure3_cfg["methods"]

    _plot_hit_rates(hit_rates, group_labels, methods, top_ks, thresholds, output_path)
    _plot_timing(timings, methods, timing_output_path)

    log.info("figure3_cross_modal.replot: wrote %s and %s", output_path, timing_output_path)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-path", required=False, help="not needed with --replot")
    ap.add_argument("--links", required=False, help="not needed with --replot")
    ap.add_argument("--config", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--timing-output", required=True)
    ap.add_argument(
        "--replot", action="store_true",
        help="redraw from --output/--timing-output's existing .json files instead of recomputing",
    )
    args = ap.parse_args()

    import yaml

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    if args.replot:
        replot(
            hit_rates_json_path=str(Path(args.output).with_suffix(".json")),
            timing_json_path=str(Path(args.timing_output).with_suffix(".json")),
            figure3_cfg=cfg["figure3"],
            mibig_highlight_classes=cfg["mibig_highlight_classes"],
            output_path=args.output,
            timing_output_path=args.timing_output,
        )
        return

    if not args.db_path or not args.links:
        ap.error("--db-path and --links are required unless --replot is given")

    run(
        db_path=args.db_path,
        links_path=args.links,
        figure3_cfg=cfg["figure3"],
        mibig_highlight_classes=cfg["mibig_highlight_classes"],
        retrieval_cfg=cfg["retrieval"],
        tanimoto_cfg=cfg["tanimoto"],
        output_path=args.output,
        timing_output_path=args.timing_output,
    )


if __name__ == "__main__":
    main()
