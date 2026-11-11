"""Figure 4: BGC-compound pair score distributions vs. random pairs.

For every MIBiG BGC with a known product compound (best-scoring compound per BGC,
same "pick one compound per BGC" rule as figures 1/2), computes two scores against
its ground truth: the database's own fingerprint cosine similarity, and the
normalized pairwise-alignment score. Both are also computed for an equal-sized set of
random (non-matching) BGC-compound pairs, as a null baseline. Real pairs should be
clearly shifted higher than random ones if RetroMol's representation is doing its job.
"""

import argparse
import json
import logging
import random
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
from tqdm import tqdm

from common import build_discovery_context, figure_style, normalize_for_alignment, open_db
from retrieval_engine import build_eval_bgcs, self_alignment_score, score_candidate_alignment, CompoundRecord


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    denom = float(np.linalg.norm(a) * np.linalg.norm(b))
    return float(np.dot(a, b) / denom) if denom > 0 else 0.0


def run(
    db_path: str,
    links_path: str,
    retrieval_cfg: dict,
    n_random_pairs: int | None,
    random_seed: int,
    output_path: str,
) -> None:
    figure_style()

    with open(links_path) as fh:
        links = json.load(fh)

    ctx = build_discovery_context()
    db = open_db(db_path)
    try:
        groups = build_eval_bgcs(db, links, {})
        queries = groups["all"]

        # Ground-truth compound fingerprints/primary_sequences, batch-fetched once.
        all_gt_ids = sorted({gt for q in queries for gt in q.ground_truth_ids})
        gt_entries = {e.id: e for e in db.get_entries(all_gt_ids)}

        real_cosine: list[float] = []
        real_alignment: list[float] = []
        for bq in tqdm(queries, desc="figure4[real pairs]", unit="bgc"):
            query_seq_norm = [normalize_for_alignment(n, ctx) for n in bq.primary_sequence]
            self_score = self_alignment_score(query_seq_norm, ctx)

            best_cos, best_align = None, None
            for gt_id in bq.ground_truth_ids:
                entry = gt_entries.get(gt_id)
                if entry is None:
                    continue
                cos = _cosine(bq.fingerprint, np.array(entry.fingerprint, dtype=np.float32))
                best_cos = cos if best_cos is None else max(best_cos, cos)

                cand = CompoundRecord(entry.id, entry.raw, entry.primary_sequence)
                pct = score_candidate_alignment(ctx, query_seq_norm, self_score, retrieval_cfg["score_mode"], cand)
                if pct is not None:
                    best_align = pct if best_align is None else max(best_align, pct)

            if best_cos is not None:
                real_cosine.append(best_cos)
            if best_align is not None:
                real_alignment.append(best_align)

        # Random (non-matching) pairs, sampled from the full compound universe.
        n_random = n_random_pairs if n_random_pairs is not None else len(queries)
        rng = random.Random(random_seed)
        all_compound_ids = [
            r[0] for r in db.con.execute("SELECT id FROM entries WHERE type = 'compound'").fetchall()
        ]

        ground_truth_by_bgc = {q.entry_id: set(q.ground_truth_ids) for q in queries}
        random_cosine: list[float] = []
        random_alignment: list[float] = []

        picks = 0
        attempts = 0
        pbar = tqdm(total=n_random, desc="figure4[random pairs]", unit="pair")
        while picks < n_random and attempts < n_random * 20:
            attempts += 1
            bq = rng.choice(queries)
            compound_id = rng.choice(all_compound_ids)
            if compound_id in ground_truth_by_bgc.get(bq.entry_id, set()):
                continue

            entry = db.get_entry(compound_id)
            if entry is None:
                continue

            query_seq_norm = [normalize_for_alignment(n, ctx) for n in bq.primary_sequence]
            self_score = self_alignment_score(query_seq_norm, ctx)

            cos = _cosine(bq.fingerprint, np.array(entry.fingerprint, dtype=np.float32))
            random_cosine.append(cos)

            cand = CompoundRecord(entry.id, entry.raw, entry.primary_sequence)
            pct = score_candidate_alignment(ctx, query_seq_norm, self_score, retrieval_cfg["score_mode"], cand)
            if pct is not None:
                random_alignment.append(pct)

            picks += 1
            pbar.update(1)
        pbar.close()
    finally:
        db.close()

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    for ax, real_vals, random_vals, title, xlabel in [
        (axes[0], real_cosine, random_cosine, "Fingerprint cosine similarity", "Cosine similarity"),
        (axes[1], real_alignment, random_alignment, "Normalized alignment score", "Normalized alignment score"),
    ]:
        x = np.linspace(0.0, 1.0, 400)
        for vals, label, color in [(real_vals, "Real BGC-compound pairs", "#4C72B0"), (random_vals, "Random pairs", "#808080")]:
            if len(vals) < 2 or np.std(vals) == 0:
                continue
            density = gaussian_kde(vals)(x)
            ax.plot(x, density, color=color, label=f"{label} (n={len(vals)})")
            ax.fill_between(x, density, color=color, alpha=0.12)
        ax.set_xlim(0, 1)
        ax.set_title(title)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Density")
        ax.legend(loc="upper left", fontsize=8)

    fig.suptitle("Figure 4 - BGC-compound pair score distributions vs. random pairs")
    fig.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)

    with open(output_path.with_suffix(".json"), "w") as fh:
        json.dump(
            {
                "real_cosine": real_cosine, "random_cosine": random_cosine,
                "real_alignment": real_alignment, "random_alignment": random_alignment,
            },
            fh,
        )

    logging.getLogger(__name__).info("figure4_score_distributions: wrote %s", output_path)


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
        retrieval_cfg=cfg["retrieval"],
        n_random_pairs=cfg["figure4"]["n_random_pairs"],
        random_seed=cfg["figure4"]["random_pairs_seed"],
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
