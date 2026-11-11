"""FASTA read/write helpers (vendored from parasect.core.parsing / parasect.core.writers)."""

from pathlib import Path


def parse_fasta_file(path_in: str | Path) -> dict[str, str]:
    """ Read a fasta file into a {header: sequence} dict. """
    path_in = Path(path_in)
    if not path_in.exists():
        raise FileNotFoundError(f"file not found: {path_in}")

    fasta: dict[str, str] = {}
    header = ""
    sequence: list[str] = []

    with open(path_in, "r") as fo:
        for line in fo:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    fasta[header] = "".join(sequence)
                header = line[1:]
                sequence = []
            else:
                sequence.append(line)

        if header:
            fasta[header] = "".join(sequence)

    return fasta


def write_fasta_file(fasta: dict[str, str], path_out: str | Path) -> None:
    """ Write a {header: sequence} dict to a fasta file, headers sorted for determinism. """
    with open(path_out, "w") as fo:
        for header in sorted(fasta):
            fo.write(f">{header}\n{fasta[header]}\n")
