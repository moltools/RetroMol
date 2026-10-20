"""MUSCLE (v3) command-line wrapper (vendored from parasect.core.muscle).

Requires the `muscle` binary on PATH -- specifically MUSCLE v3's
`-in1/-in2/-profile` interface, which MUSCLE v5 does not support.
"""

import subprocess
from pathlib import Path


DEFAULT_TIMEOUT_S = 180


def run_muscle_profile(
    path_in: str | Path, path_in_alignment: str | Path, path_out: str | Path, timeout: float = DEFAULT_TIMEOUT_S
) -> None:
    """
    Align a single new sequence against an existing MUSCLE alignment via profile-profile alignment.

    :param path_in: fasta file with the (single) new sequence to align.
    :param path_in_alignment: fasta file with the existing reference alignment.
    :param path_out: path to write the resulting alignment to.
    :param timeout: seconds to wait before giving up -- a stuck `muscle` process (seen under
        Rosetta emulation of the old 3.8.1551 build on Apple Silicon) would otherwise hang
        forever instead of raising.
    :raises subprocess.TimeoutExpired: if `muscle` doesn't finish within `timeout` seconds.
    """
    command = [
        "muscle", "-quiet", "-profile",
        "-in1", str(path_in_alignment),
        "-in2", str(path_in),
        "-out", str(path_out),
    ]
    subprocess.run(command, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, timeout=timeout)
