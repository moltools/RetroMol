#!/usr/bin/env python3

import numpy as np
import pandas as pd

from retromol_alignment.scoring import create_scoring_matrix
from retromol_alignment.aligner import setup_aligner
from retromol_alignment.pairwise import Converter, align
from retromol_alignment.msa import calculate_msa


def main() -> None:
    # Dummy scoring matrix with three items: A, B, C
    df = pd.DataFrame(
        {
            "A": [1, -1, -2],
            "B": [-1, 1, -1],
            "C": [-2, -1, 1]
        },
        index=["A", "B", "C"],
    )

    m = create_scoring_matrix(df)

    aligner = setup_aligner(m, mode="global")

    seq_a = ["A", "B", "C", "A"]
    seq_b = ["A", "C", "B", "A"]
    seq_c = ["A", "B", "A"]

    converter = Converter(
        to_identifier=lambda x: m.alphabet.index(x),
        from_identifier=lambda i: m.alphabet[i],
        gap_repr=np.int32(-1),
        insert_repr=np.int32(-2)
    )

    seq_a_arr = converter.to_int_array(seq_a)
    print(seq_a_arr)

    score, ala, alb = align(aligner, seq_a, seq_b, converter)
    print(score)
    print(ala)
    print(alb)

    msa, new_order = calculate_msa(aligner, [seq_a, seq_b, seq_c], converter, 0, [True, True, True])
    names = [f"seq_{i+1}" for i in new_order]
    width = 1
    gap_repr = ' '
    d = {x: x for x in m.alphabet}
    for i, (s, is_rev, al) in enumerate(msa):
        s = f"{names[i]} ({s:>4.2f}; {'rev' if is_rev else 'fwd'}) > " + "|".join(f"{d.get(x, gap_repr)!s:>{width}}" for x in al)
        print(s)


if __name__ == "__main__":
    main()
