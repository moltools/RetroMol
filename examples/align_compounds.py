# -*- coding: utf-8 -*-

"""Example script that shows how to align primary sequences of two compounds using RetroMol."""

from typing import List, Tuple

from versalign.aligner import setup_aligner
from versalign.msa import calc_msa
from versalign.printing import format_alignment
from versalign.scoring import create_substituion_matrix_dynamically

from retromol.api import run_retromol_with_timeout
from retromol.io import Input as RetroMolInput
from retromol.readout import linear_readout


def main() -> None:
    """Main function to align two compounds."""
    
    # First parse both compounds with RetroMol
    inp1 = RetroMolInput("dictyostatin", r"C[C@H]1CC[C@H]([C@@H]([C@@H](OC(=O)/C=C\C=C\[C@H]([C@H](C[C@@H](/C=C\[C@@H]([C@@H]([C@H](C1)C)O)C)O)O)C)[C@@H](C)/C=C\C=C)C)O")
    inp2 = RetroMolInput("discodermolide", r"C[C@H]1[C@@H](OC(=O)[C@@H]([C@H]1O)C)C[C@@H](/C=C\[C@H](C)[C@@H]([C@@H](C)/C=C(/C)\C[C@H](C)[C@H]([C@H](C)[C@H]([C@@H](C)/C=C\C=C)OC(=O)N)O)O)O")
    res1 = run_retromol_with_timeout(inp1)  # uses default ruleset
    res2 = run_retromol_with_timeout(inp2)  # uses default ruleset
    cov1 = res1.best_total_coverage()
    cov2 = res2.best_total_coverage()
    print(f"Coverage for dictyostatin: {cov1:.1%}")
    print(f"Coverage for discodermolide: {cov2:.1%}")

    # Get linear readouts from both compounds
    readout1 = linear_readout(res1, require_identified=True)
    readout2 = linear_readout(res2, require_identified=True)
    readouts = [("dictyostatin", readout1), ("discodermolide", readout2)]

    # Extract primary sequences from readouts
    records: List[Tuple[str, str]] = []
    for label, readout in readouts:
        for level_idx, level in enumerate(readout["levels"]):
            for path_idx, path in enumerate(level["strict_paths"]):
                path = path["ordered_monomers"]
                if len(path) <= 3: continue  # skip too short paths
                seq_fwd = [m["identity"] for m in path]
                seq_rev = list(reversed(seq_fwd))
                records.append((f"{label}_{level_idx}_{path_idx}_fwd", seq_fwd))
                records.append((f"{label}_{level_idx}_{path_idx}_rev", seq_rev))

    # Unzip labels and sequences
    labels, seqs = zip(*records)

    # Align sequences using versalign
    objs = list(set([x for seq in seqs for x in seq])) + ["-"]  # include gap character
    objs.sort()
    sm, _ = create_substituion_matrix_dynamically(objs)
    aligner = setup_aligner(sm, mode="global")
    msa, order = calc_msa(aligner, seqs, gap_repr="-")
    reordered_labels = [labels[i] for i in order]
    print("\n" + format_alignment(msa, names=reordered_labels))

    # NOTE: see versalign documentation for more alignment options, including building 
    #       more elaborate substitution matrices.


if __name__ == "__main__":
    main()
