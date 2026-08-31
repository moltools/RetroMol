"""Figure 1: RetroMol coverage distributions.

Left panel: all NPAtlas compounds vs. NPAtlas compounds restricted to a handful of
NPClassifier chemical classes (default: oligopeptides, macrolides, cyclic
polyketides), meant to show coverage is especially good for those classes.

Right panel: all MIBiG BGCs (one point per BGC accession, the best-coverage
compound among that accession's associated product compounds, since a BGC can have
more than one) vs. the same, restricted to MIBiG's own Polyketide and Nonribosomal
peptide biosynthetic_class labels.

Both panels overlay density curves (gaussian KDE) on the same [0, 1] coverage axis
so they're directly comparable, plus the raw sample counts in the legend.

Coverage is read straight off `entries.coverage` (see
src/retromol_database/duckdb.py / database/scripts/load_compounds.py); no RetroMol
re-parsing needed here.
"""

import argparse
import json
import logging
from pathlib import Path

import duckdb
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde

from common import exact_label_group_accessions, figure_style

log = logging.getLogger(__name__)

_COLORS = ["#4C72B0", "#DD8452", "#55A868", "#C44E52"]


def _kde_line(values: list[float], ax, label: str, color: str, linestyle: str = "-") -> None:
    x = np.linspace(0.0, 1.0, 400)
    if len(values) < 2 or np.std(values) == 0:
        ax.axvline(values[0] if values else 0.0, color=color, linestyle=linestyle, label=f"{label} (n={len(values)})")
        return
    density = gaussian_kde(values)(x)
    ax.plot(x, density, color=color, linestyle=linestyle, label=f"{label} (n={len(values)})")
    ax.fill_between(x, density, color=color, alpha=0.12)


def _npatlas_panel(ax, con, coverage_by_entry: dict[str, float], highlight_classes: dict) -> None:
    npatlas_ids = {
        r[0]
        for r in con.execute(
            "SELECT DISTINCT entry_id FROM entry_sources WHERE database_name = 'NPAtlas'"
        ).fetchall()
    }
    npatlas_entries = {eid: cov for eid, cov in coverage_by_entry.items() if eid in npatlas_ids}
    all_values = list(npatlas_entries.values())
    _kde_line(all_values, ax, "All NPAtlas compounds", "#808080", linestyle="--")

    for i, (_, group) in enumerate(highlight_classes.items()):
        labels = group["npclassifier_labels"]
        level = group.get("npclassifier_level", "class")
        placeholders = ",".join(["?"] * len(labels))
        rows = con.execute(
            f"""
            SELECT DISTINCT entry_id FROM chemical_class_annotations
            WHERE level = ? AND label IN ({placeholders})
            """,
            [level, *labels],
        ).fetchall()
        group_ids = {r[0] for r in rows} & npatlas_entries.keys()
        values = [npatlas_entries[eid] for eid in group_ids]
        _kde_line(values, ax, group["label"], _COLORS[i % len(_COLORS)])

    ax.set_title("NPAtlas compounds")
    ax.set_xlabel("RetroMol coverage")
    ax.set_ylabel("Density")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper left", fontsize=8)


def _mibig_panel(
    ax,
    con,
    coverage_by_entry: dict[str, float],
    links: dict[str, dict],
    highlight_classes: dict,
) -> None:
    # One value per accession: the best-coverage compound among its associated
    # products (a BGC entry can have more than one associated product compound).
    best_by_accession: dict[str, float] = {}
    for accession, ids in links.items():
        if not ids["bgc_entry_ids"] or not ids["compound_entry_ids"]:
            continue
        covs = [coverage_by_entry[eid] for eid in ids["compound_entry_ids"] if eid in coverage_by_entry]
        if covs:
            best_by_accession[accession] = max(covs)

    all_values = list(best_by_accession.values())
    _kde_line(all_values, ax, "All MIBiG BGCs", "#808080", linestyle="--")

    group_label_sets = {key: frozenset(g["biosynthetic_class_labels"]) for key, g in highlight_classes.items()}
    group_accessions_by_key = exact_label_group_accessions(con, links, group_label_sets)

    for i, (key, group) in enumerate(highlight_classes.items()):
        values = [best_by_accession[a] for a in group_accessions_by_key[key] if a in best_by_accession]
        _kde_line(values, ax, group["label"], _COLORS[i % len(_COLORS)])

    ax.set_title("MIBiG BGCs (best-scoring product per BGC)")
    ax.set_xlabel("RetroMol coverage")
    ax.set_xlim(0, 1)
    ax.legend(loc="upper left", fontsize=8)


def run(
    db_path: str,
    links_path: str,
    npatlas_highlight_classes: dict,
    mibig_highlight_classes: dict,
    output_path: str,
) -> None:
    figure_style()

    with open(links_path) as fh:
        links = json.load(fh)

    con = duckdb.connect(db_path, read_only=True)
    try:
        coverage_by_entry = {
            r[0]: r[1]
            for r in con.execute(
                "SELECT id, coverage FROM entries WHERE type = 'compound' AND coverage IS NOT NULL"
            ).fetchall()
        }

        fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
        _npatlas_panel(axes[0], con, coverage_by_entry, npatlas_highlight_classes)
        _mibig_panel(axes[1], con, coverage_by_entry, links, mibig_highlight_classes)
    finally:
        con.close()

    fig.suptitle("Figure 1 - RetroMol coverage distributions")
    fig.tight_layout()

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)
    log.info("figure1_coverage: wrote %s", output_path)


def main() -> None:
    logging.basicConfig(level=logging.INFO)
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db-path", required=True)
    ap.add_argument("--links", required=True)
    ap.add_argument("--config", required=True, help="path to benchmark/config.yaml")
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    import yaml

    with open(args.config) as fh:
        cfg = yaml.safe_load(fh)

    run(
        db_path=args.db_path,
        links_path=args.links,
        npatlas_highlight_classes=cfg["npatlas_highlight_classes"],
        mibig_highlight_classes=cfg["mibig_highlight_classes"],
        output_path=args.output,
    )


if __name__ == "__main__":
    main()
