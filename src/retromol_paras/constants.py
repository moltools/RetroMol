"""Data-file paths and parsed constants for retromol_paras (vendored from parasect.core.constants)."""

from importlib.resources import files
from pathlib import Path

from retromol_paras.tabular import Tabular

DATA_DIR = Path(files("retromol_paras.data").joinpath(""))


def _get_path(file_name: str) -> Path:
    return DATA_DIR / file_name


def _read_positions(path_in: Path, start_position: int) -> list[int]:
    """ Parse a tab-separated list of reference positions, offset by `start_position`. """
    text = path_in.read_text().strip()
    return [int(p) - start_position for p in text.split("\t")]


def _parse_amino_acid_properties(path_in: Path) -> dict[str, list[float]]:
    """ Parse the 15-dimensional physicochemical property vector per amino acid (+ gap, "-"). """
    data = Tabular(path_in, separator="\t")
    properties: dict[str, list[float]] = {}
    for row_id in data.rows:
        values = [float(v) for v in data.get_row_values(row_id)[1:]]
        if len(values) != 15:
            raise ValueError(f"amino acid property vector for {row_id} is not of length 15")
        properties[data.get_row_value(row_id, "AA")] = values
    return properties


ALIGNMENT_FILE = _get_path("structure_alignment.fasta")
ALIGNMENT_REFERENCE_ID = "BAA00406.1.A1"

HMM2_FILE = _get_path("AMP-binding_hmmer2.hmm")
HMM3_FILE = _get_path("AMP-binding_full.hmm")

# HMMER2-hit-relative signature positions (offset 0: parasect's stachelhaus_hmm2.txt/
# active_site_hmm2.txt are already given relative to the HMM2 hit start).
HMM2_POSITIONS_SIGNATURE = _read_positions(_get_path("stachelhaus_hmm2.txt"), 0)
HMM2_POSITIONS_EXTENDED_SIGNATURE = _read_positions(_get_path("active_site_hmm2.txt"), 0)
HMM2_POSITION_LYSINE = [36]

# Profile-alignment (MUSCLE, reference-sequence-relative) signature positions -- offset 66
# because parasect's stachelhaus.txt/active_site.txt give absolute 1AMU reference positions.
POSITIONS_SIGNATURE = _read_positions(_get_path("stachelhaus.txt"), 66)
POSITIONS_EXTENDED_SIGNATURE = _read_positions(_get_path("active_site.txt"), 66)

PROPERTIES = _parse_amino_acid_properties(_get_path("physicochemical_properties.txt"))

SMILES_FILE = _get_path("database_files/smiles.tsv")
PARASECT_DATASET_FILE = _get_path("database_files/parasect_dataset.txt")
INCLUDED_SUBSTRATES_FILE = _get_path("included_substrates.txt")
INCLUDED_SUBSTRATES_FILE_BACTERIAL = _get_path("included_substrates_bacterial.txt")
