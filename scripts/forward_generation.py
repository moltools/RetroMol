#!/usr/bin/env python3
from retromol.forward import forward_generation

def main() -> None:
    motifs = ["ethanoic acid", "A1", "B1", "alanine", "C1", "C1", "C1", "D1", "glycine", "proline", "tryptophan"]
    smiles = forward_generation(motifs)
    print(smiles)

if __name__ == "__main__":
    main()