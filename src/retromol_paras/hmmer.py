"""HMMER (v2 + v3) command-line wrappers and result parsing (vendored from parasect.core.hmmer).

Requires `hmmpfam2` (HMMER2) and `hmmscan` + `hmmpress` (HMMER3) on PATH.
"""

import fcntl
import subprocess
from pathlib import Path

from Bio import SearchIO
from Bio.SearchIO._model import HSP

from retromol_paras.fasta import parse_fasta_file

DEFAULT_TIMEOUT_S = 180


def ensure_hmm3_pressed(hmm3_file: str | Path, timeout: float = DEFAULT_TIMEOUT_S) -> None:
    """
    Run `hmmpress` on the packaged HMM3 profile once, if its binary index files aren't
    there yet. Locked (same pattern as retromol_paras.train's training lock) since this
    file lives at a single shared package-data path -- concurrent worker processes each
    seeing "not pressed yet" on their first call and racing to press it simultaneously
    would corrupt the index files, not just waste work.

    :raises subprocess.TimeoutExpired: if `hmmpress` doesn't finish within `timeout` seconds.
    """
    hmm3_file = Path(hmm3_file)
    if all((hmm3_file.parent / f"{hmm3_file.name}.{ext}").exists() for ext in ("h3f", "h3i", "h3m", "h3p")):
        return

    with open(hmm3_file.parent / f"{hmm3_file.name}.press.lock", "w") as lock_fo:
        fcntl.flock(lock_fo.fileno(), fcntl.LOCK_EX)
        try:
            if all((hmm3_file.parent / f"{hmm3_file.name}.{ext}").exists() for ext in ("h3f", "h3i", "h3m", "h3p")):
                return  # another process pressed it while we waited for the lock
            subprocess.run(
                ["hmmpress", str(hmm3_file)], check=True,
                stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout,
            )
        finally:
            fcntl.flock(lock_fo.fileno(), fcntl.LOCK_UN)


def run_hmmscan(hmm_file: str | Path, fasta_file: str | Path, out_file: str | Path, timeout: float = DEFAULT_TIMEOUT_S) -> None:
    """ Run HMMER3's `hmmscan` on a fasta file against a (pressed) HMM database.

    :raises subprocess.TimeoutExpired: if `hmmscan` doesn't finish within `timeout` seconds.
    """
    with open(out_file, "w") as out:
        subprocess.run(["hmmscan", str(hmm_file), str(fasta_file)], check=True, stdout=out, timeout=timeout)


def run_hmmpfam2(hmm_file: str | Path, fasta_file: str | Path, out_file: str | Path, timeout: float = DEFAULT_TIMEOUT_S) -> None:
    """ Run HMMER2's `hmmpfam2` on a fasta file against an HMM database.

    :raises subprocess.TimeoutExpired: if `hmmpfam2` doesn't finish within `timeout` seconds.
    """
    with open(out_file, "w") as out:
        subprocess.run(["hmmpfam2", str(hmm_file), str(fasta_file)], check=True, stdout=out, timeout=timeout)


def parse_hmm_results(path_in: str | Path, hmmer_version: int) -> dict[str, HSP]:
    """
    Parse hmmscan/hmmpfam2 output, keeping AMP-binding / AMP-binding_C hits with bitscore > 20.

    :param path_in: path to the hmmscan (v3) or hmmpfam2 (v2) output file.
    :param hmmer_version: 2 or 3.
    :return: {"<query_id>|<hit_id>|<query_start>-<query_end>": HSP}
    """
    if hmmer_version not in (2, 3):
        raise ValueError(f"unknown HMMer version: {hmmer_version}")

    filtered_hits: dict[str, HSP] = {}
    hmmer_string = f"hmmer{hmmer_version}-text"

    for result in SearchIO.parse(str(path_in), hmmer_string):
        for hsp in result.hsps:
            if hsp.bitscore > 20 and hsp.hit_id in ("AMP-binding", "AMP-binding_C"):
                header = f"{result.id}|{hsp.hit_id}|{hsp.query_start}-{hsp.query_end}"
                filtered_hits[header] = hsp

    return filtered_hits


def rename_sequences(path_in: str | Path, path_out_dir: str | Path) -> tuple[Path, Path]:
    """
    Rename fasta headers to plain integers before running HMMER (some headers break HMMER's parser),
    returning (mapping_file, renamed_fasta_file). Reverse with `reverse_renaming`.
    """
    path_in = Path(path_in)
    if not path_in.exists():
        raise FileNotFoundError(f"file not found: {path_in}")

    path_out_dir = Path(path_out_dir)
    path_out_dir.mkdir(parents=True, exist_ok=True)

    mapping_file = path_out_dir / "mapping.txt"
    renamed_fasta_file = path_out_dir / "renamed_fasta.fasta"

    id_to_seq = parse_fasta_file(path_in)

    with open(renamed_fasta_file, "w") as fo_fasta, open(mapping_file, "w") as fo_mapping:
        for counter, (seq_id, seq) in enumerate(id_to_seq.items(), start=1):
            fo_fasta.write(f">{counter}\n{seq}\n")
            fo_mapping.write(f"{counter}\t{seq_id}\n")

    return mapping_file, renamed_fasta_file


def load_renaming_map(mapping_file: str | Path) -> dict[str, str]:
    """ Load the renamed-id -> original-id mapping written by `rename_sequences`. """
    new_to_original: dict[str, str] = {}
    with open(mapping_file, "r") as fo:
        for line in fo:
            new, _, original = line.rstrip("\n").partition("\t")
            new_to_original[new] = original
    return new_to_original
