"""
TODO:
- parse compounds and make sure coverage is high
- parse gbks; get distribution of classes
- embed true pairs vs. random pairs and compare distribution of cosine scores
- align true paris vs. random pairs and compare distribution of alignment scores
"""

import argparse
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from rdkit import RDLogger
from tqdm import tqdm

from retromol.model.submission import Submission
from retromol.model.rules import RuleSet
from retromol.pipelines.parsing import run_retromol_with_timeout


RDLogger.DisableLog("rdApp.*")


def cli() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--work-dir", type=str, required=True)
    parser.add_argument("--compounds-xlsx", required=True)
    return parser.parse_args()


def main() -> None:
    args = cli()

    work_dir = Path(args.work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    compounds = Path(args.compounds_xlsx)
    compounds_df = pd.read_excel(compounds)
    unique_pairs = compounds_df[["Cluster", "True SMILES"]].drop_duplicates()

    rules = RuleSet.load_default(match_stereochemistry=False)

    # Parse compounds with RetroMol
    coverages = []
    failed_smis = []

    out_fn = work_dir / "prism4_gold_standard.tsv"
    with open(out_fn, "w") as fo:
        fo.write("smiles\tcoverage\tfailed\n")

        for row in tqdm(unique_pairs.iterrows()):
            cluster_fn, smiles = row[1]
            submission = Submission(name=cluster_fn, smiles=smiles)
            try:
                failed = False
                result = run_retromol_with_timeout(submission=submission, rules=rules)
                coverage = result.calculate_coverage()
                coverages.append(coverage)
            except:
                failed = True
                failed_smis.append(smiles)
                coverage = 0.0
                coverages.append(coverage)

            fo.write(f"{smiles}\t{coverage}\t{failed}\n")
            fo.flush()

    # Plot coverages in barplot (20 bins; range 0-1)
    plt.hist(coverages, bins=20, range=(0, 1))
    plt.xlabel("Coverage")
    plt.ylabel("Frequency")
    plt.title("Distribution of RetroMol Coverage for PRISM4 Gold Standard")
    plt.savefig(work_dir / "coverage_distribution.png", dpi=300)
    plt.close()

    for smiles in failed_smis:
        print(f"Failed to parse: {smiles}")
    print(f"Total failed: {len(failed_smis)}")


if __name__ == "__main__":
    main()
