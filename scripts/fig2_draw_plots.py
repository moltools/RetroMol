import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tsv", required=True, type=Path)
    parser.add_argument("--outdir", default="plots", type=Path)

    parser.add_argument(
        "--polyketide-classes",
        type=Path,
        default=None,
        help="Optional txt file with one np_superclass per line to define polyketide compounds.",
    )
    parser.add_argument(
        "--peptide-classes",
        type=Path,
        default=None,
        help="Optional txt file with one np_superclass per line to define peptide compounds.",
    )

    return parser.parse_args()


def read_class_file(path: Path) -> set[str]:
    classes = set()

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                classes.add(line.lower())

    return classes


def split_superclasses(value: str) -> list[str]:
    return [x.strip() for x in str(value).split(";") if x.strip()]


def make_mask_from_class_file(df: pd.DataFrame, class_file: Path) -> pd.Series:
    allowed_classes = read_class_file(class_file)

    return df["np_superclass"].apply(
        lambda value: any(
            superclass.lower() in allowed_classes
            for superclass in split_superclasses(value)
        )
    )


def coverage_threshold_counts(values: pd.Series) -> dict[str, int]:
    values = pd.to_numeric(values, errors="coerce").dropna()
    values = values[(values >= 0) & (values <= 1)]

    return {
        "n": len(values),
        "le_0.1": int((values <= 0.1).sum()),
        "ge_0.5": int((values >= 0.5).sum()),
        "ge_0.9": int((values >= 0.9).sum()),
    }


def make_histogram(values: pd.Series, title: str, out_file: Path) -> dict[str, int]:
    values = pd.to_numeric(values, errors="coerce").dropna()
    values = values[(values >= 0) & (values <= 1)]

    counts = coverage_threshold_counts(values)

    bins = np.linspace(0, 1, 21)

    plt.figure(figsize=(6, 4.8))
    plt.hist(values, bins=bins, color="#56b4e9")
    plt.xlim(0, 1)
    plt.xlabel("Coverage")
    plt.ylabel("Count")
    plt.title(f"{title} (n={counts['n']})")

    threshold_text = (
        f"coverage <= 0.1: {counts['le_0.1']} compounds\n"
        f"coverage >= 0.5: {counts['ge_0.5']} compounds\n"
        f"coverage >= 0.9: {counts['ge_0.9']} compounds"
    )

    plt.figtext(
        0.5,
        0.01,
        threshold_text,
        ha="center",
        va="bottom",
        fontsize=9,
    )

    plt.tight_layout(rect=(0, 0.14, 1, 1))
    plt.savefig(out_file, dpi=300)
    plt.close()

    return counts


def superclass_coverage_summary(df: pd.DataFrame, mask: pd.Series | None = None) -> pd.DataFrame:
    if mask is None:
        sub = df[["np_superclass", "coverage"]].copy()
    else:
        sub = df.loc[mask, ["np_superclass", "coverage"]].copy()

    exploded = (
        sub.dropna(subset=["np_superclass"])
        .assign(np_superclass=lambda x: x["np_superclass"].astype(str).str.split(";"))
        .explode("np_superclass")
    )

    exploded["np_superclass"] = exploded["np_superclass"].str.strip()
    exploded = exploded[exploded["np_superclass"] != ""]
    exploded["coverage"] = pd.to_numeric(exploded["coverage"], errors="coerce")
    exploded = exploded.dropna(subset=["coverage"])
    exploded = exploded[(exploded["coverage"] >= 0) & (exploded["coverage"] <= 1)]

    return (
        exploded.groupby("np_superclass")["coverage"]
        .agg(
            count="count",
            mean_coverage="mean",
            std_coverage="std",
        )
        .reset_index()
    )


def write_superclass_summary(summary: pd.DataFrame, out_file: Path, sort_by: str) -> None:
    summary = summary.sort_values(sort_by, ascending=False)

    with open(out_file, "w") as f:
        f.write("np_superclass\tcount\tmean_coverage\tstd_coverage\n")

        for _, row in summary.iterrows():
            std = row["std_coverage"]
            std_str = "" if pd.isna(std) else f"{std:.6f}"

            f.write(
                f"{row['np_superclass']}\t"
                f"{int(row['count'])}\t"
                f"{row['mean_coverage']:.6f}\t"
                f"{std_str}\n"
            )


def write_plot_summary(summary: dict[str, dict[str, int]], out_file: Path) -> None:
    with open(out_file, "w") as f:
        f.write("plot\tn\tcoverage_le_0.1\tcoverage_ge_0.5\tcoverage_ge_0.9\n")

        for plot_name, counts in summary.items():
            f.write(
                f"{plot_name}\t"
                f"{counts['n']}\t"
                f"{counts['le_0.1']}\t"
                f"{counts['ge_0.5']}\t"
                f"{counts['ge_0.9']}\n"
            )


def main() -> None:
    args = cli()
    args.outdir.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.tsv, sep="\t")

    if "coverage" not in df.columns:
        raise ValueError("Input TSV must contain a 'coverage' column")

    if "np_superclass" not in df.columns:
        raise ValueError("Input TSV must contain an 'np_superclass' column")

    unique_np_superclass = sorted([str(x) for x in df["np_superclass"].unique().tolist()])
    for superclass in unique_np_superclass:
        print(superclass)

    df["coverage"] = pd.to_numeric(df["coverage"], errors="coerce")
    df["np_superclass"] = df["np_superclass"].fillna("").astype(str)

    if args.polyketide_classes is not None:
        polyketide_mask = make_mask_from_class_file(df, args.polyketide_classes)
    else:
        polyketide_mask = df["np_superclass"].str.contains(
            "polyketide",
            case=False,
            regex=False,
        )

    if args.peptide_classes is not None:
        peptide_mask = make_mask_from_class_file(df, args.peptide_classes)
    else:
        peptide_mask = df["np_superclass"].str.contains(
            "peptide",
            case=False,
            regex=False,
        )

    plot_summary = {}

    plot_summary["all"] = make_histogram(
        df["coverage"],
        title="Coverage distribution - all compounds",
        out_file=args.outdir / "coverage_all.svg",
    )

    plot_summary["polyketide"] = make_histogram(
        df.loc[polyketide_mask, "coverage"],
        title="Coverage distribution - polyketide compounds",
        out_file=args.outdir / "coverage_polyketide.svg",
    )

    plot_summary["peptide"] = make_histogram(
        df.loc[peptide_mask, "coverage"],
        title="Coverage distribution - peptide compounds",
        out_file=args.outdir / "coverage_peptide.svg",
    )

    polyketide_summary = superclass_coverage_summary(df, polyketide_mask)
    peptide_summary = superclass_coverage_summary(df, peptide_mask)
    all_summary = superclass_coverage_summary(df)

    write_superclass_summary(
        polyketide_summary,
        args.outdir / "np_superclasses_polyketide.txt",
        sort_by="count",
    )

    write_superclass_summary(
        peptide_summary,
        args.outdir / "np_superclasses_peptide.txt",
        sort_by="count",
    )

    write_superclass_summary(
        all_summary,
        args.outdir / "np_superclasses_by_mean_coverage.txt",
        sort_by="mean_coverage",
    )

    write_plot_summary(
        plot_summary,
        args.outdir / "coverage_plot_summary.txt",
    )

    for plot_name, counts in plot_summary.items():
        print(
            f"{plot_name}: "
            f"n={counts['n']}, "
            f"coverage <= 0.1: {counts['le_0.1']}, "
            f"coverage >= 0.5: {counts['ge_0.5']}, "
            f"coverage >= 0.9: {counts['ge_0.9']}"
        )

    print(f"Wrote plots and summaries to: {args.outdir}")


if __name__ == "__main__":
    main()