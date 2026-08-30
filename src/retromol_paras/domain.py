"""Adenylation-domain signature extraction (vendored from parasect.core.domain)."""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from Bio.SearchIO._model.hsp import HSP

from retromol_paras.constants import (
    ALIGNMENT_FILE,
    ALIGNMENT_REFERENCE_ID,
    HMM2_POSITION_LYSINE,
    HMM2_POSITIONS_EXTENDED_SIGNATURE,
    HMM2_POSITIONS_SIGNATURE,
    POSITIONS_EXTENDED_SIGNATURE,
    POSITIONS_SIGNATURE,
)
from retromol_paras.fasta import parse_fasta_file
from retromol_paras.muscle import run_muscle_profile

VALID_RESIDUES = set("ACDEFGHIKLMNPQRSTVWY-")


def _reference_positions_in_alignment(reference_positions: list[int], aligned_reference: str) -> list[int]:
    """ Map ungapped reference positions onto their (gapped) index within an aligned reference sequence. """
    positions = []
    ungapped_position = 0
    for idx, residue in enumerate(aligned_reference):
        if residue != "-":
            if ungapped_position in reference_positions:
                positions.append(idx)
            ungapped_position += 1
    return positions


def _reference_positions_in_hmm_alignment(
    query_sequence: str, reference_sequence: str, reference_positions: list[int]
) -> Optional[list[int]]:
    """
    Same as `_reference_positions_in_alignment`, but for an HMMER-produced query/profile
    alignment pair, where the "reference" is the HMM consensus track (may use "." for gaps
    too, and may be shorter than the query).
    """
    if len(reference_sequence) < len(reference_positions):
        raise ValueError(f"reference sequence too short: {len(reference_sequence)} < {len(reference_positions)}")
    if not (len(reference_sequence) <= len(query_sequence)):
        raise ValueError(f"reference sequence too long: {len(reference_sequence)} > {len(query_sequence)}")

    positions = []
    ungapped_position = 0
    for idx, residue in enumerate(reference_sequence):
        if residue in "-.":
            continue
        if ungapped_position in reference_positions:
            positions.append(idx)
        ungapped_position += 1

    if len(positions) != len(reference_positions):
        return None
    return positions


def _gap_adjusted_positions(query: str, positions: list[int], offset: int) -> list[Optional[int]]:
    """ Convert gapped-alignment column indices into ungapped sequence coordinates (offset by `offset`). """
    adjusted: list[Optional[int]] = []
    position = offset
    for i, char in enumerate(query):
        if i in positions:
            adjusted.append(position if char != "-" else None)
        if char != "-":
            position += 1
    return adjusted


def _align_to_reference(domain_name: str, domain_sequence: str, path_temp_dir: str | Path) -> tuple[str, str]:
    """ Profile-align one domain sequence against the packaged reference alignment via MUSCLE3. """
    path_temp_dir = Path(path_temp_dir)
    temp_in = path_temp_dir / "temp_in_alignment.fasta"
    temp_out = path_temp_dir / "temp_out_alignment.fasta"

    temp_in.write_text(f">{domain_name}\n{domain_sequence}\n")
    run_muscle_profile(path_in=temp_in, path_in_alignment=ALIGNMENT_FILE, path_out=temp_out)

    aligned = parse_fasta_file(temp_out)
    return aligned[domain_name], aligned[ALIGNMENT_REFERENCE_ID]


@dataclass
class AdenylationDomain:
    """ A single adenylation-domain hit, with its extracted 34-residue extended signature. """

    protein_name: str
    start: int
    end: int
    domain_nr: int = 0
    sequence: str = ""
    protein_sequence: str = ""
    signature: str = ""
    extended_signature: str = ""
    signature_positions: list[Optional[int]] = field(default_factory=list)
    extended_signature_positions: list[Optional[int]] = field(default_factory=list)

    def overlaps(self, other: "AdenylationDomain", threshold: int = 50) -> bool:
        if other.start <= self.start <= other.end:
            return other.end - self.start >= threshold
        if self.start <= other.start <= self.end:
            return self.end - other.start >= threshold
        return False

    def set_signatures_from_hmm(self, hit_n_terminal: HSP, hit_c_terminal: Optional[HSP] = None) -> None:
        """ Extract signature + extended signature from an HMMER2 N-terminal (+ optional C-terminal) hit. """
        profile = hit_n_terminal.aln[1].seq
        query = hit_n_terminal.aln[0].seq
        offset = hit_n_terminal.hit_start
        query_offset = hit_n_terminal.query_start

        signature_cols = _reference_positions_in_hmm_alignment(
            query, profile, [p - offset for p in HMM2_POSITIONS_SIGNATURE]
        )
        if signature_cols:
            signature = "".join(query[i] for i in signature_cols)
            if all(c in VALID_RESIDUES for c in signature):
                self.signature = signature
                self.signature_positions = _gap_adjusted_positions(query, signature_cols, query_offset)

        lysine = None
        if hit_c_terminal is not None:
            profile_c = hit_c_terminal.aln[1].seq
            query_c = hit_c_terminal.aln[0].seq
            offset_c = hit_c_terminal.hit_start
            lysine_cols = _reference_positions_in_hmm_alignment(
                query_c, profile_c, [p - offset_c for p in HMM2_POSITION_LYSINE]
            )
            if lysine_cols:
                lysine = query_c[lysine_cols[0]]

        if self.signature:
            self.signature += lysine if lysine and lysine in VALID_RESIDUES and lysine != "-" else "K"

        extended_cols = _reference_positions_in_hmm_alignment(
            query, profile, [p - offset for p in HMM2_POSITIONS_EXTENDED_SIGNATURE]
        )
        if extended_cols:
            extended = "".join(query[i] for i in extended_cols)
            if all(c in VALID_RESIDUES for c in extended):
                self.extended_signature = extended
                self.extended_signature_positions = _gap_adjusted_positions(query, extended_cols, query_offset)

    def set_signatures_from_profile(self, path_temp_dir: str | Path) -> None:
        """ Extract signature + extended signature via MUSCLE3 profile alignment against the reference set. """
        if not self.sequence:
            raise ValueError("domain sequence must be set before profile alignment")

        aligned_domain, aligned_reference = _align_to_reference("DOMAIN_TO_QUERY", self.sequence, path_temp_dir)

        signature_cols = _reference_positions_in_alignment(POSITIONS_SIGNATURE, aligned_reference)
        self.signature = "".join(aligned_domain[i] for i in signature_cols)
        self.signature_positions = _gap_adjusted_positions(aligned_domain, signature_cols, self.start)

        extended_cols = _reference_positions_in_alignment(POSITIONS_EXTENDED_SIGNATURE, aligned_reference)
        self.extended_signature = "".join(aligned_domain[i] for i in extended_cols)
        self.extended_signature_positions = _gap_adjusted_positions(aligned_domain, extended_cols, self.start)
