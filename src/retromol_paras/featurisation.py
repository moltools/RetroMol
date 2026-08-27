"""Extracts adenylation domains from a fasta file and featurises them (vendored from
parasect.core.featurisation), using PARAS' full default pipeline: HMMER2 as the primary
signature-extraction method, with HMMER3 (+ MUSCLE3 profile alignment) picking up any
domain HMMER2 missed.
"""

from pathlib import Path
from typing import Optional

from Bio.SearchIO._model.hsp import HSP

from retromol_paras.constants import HMM2_FILE, HMM3_FILE, PROPERTIES
from retromol_paras.domain import AdenylationDomain
from retromol_paras.fasta import parse_fasta_file
from retromol_paras.hmmer import ensure_hmm3_pressed, load_renaming_map, parse_hmm_results, rename_sequences, run_hmmpfam2, run_hmmscan

Hit = tuple[str, int, int, str]  # (hit_id, hit_start, hit_end, hit_key)


def get_domain_features(extended_signature: str) -> list[float]:
    """ Concatenate the 15-dim physicochemical property vector for each residue in an extended signature. """
    features: list[float] = []
    for residue in extended_signature:
        features.extend(PROPERTIES[residue])
    return features


def _merge_hits(hits: list[Hit]) -> Hit:
    seq_id, hit_id, _ = hits[0][3].split("|")
    for hit in hits:
        seq_id_2, hit_id_2, _ = hit[3].split("|")
        if seq_id_2 != seq_id or hit_id_2 != hit_id:
            raise ValueError("cannot merge hits from different sequences/hit types")
    hit_start = min(h[1] for h in hits)
    hit_end = max(h[2] for h in hits)
    return (hit_id, hit_start, hit_end, f"{seq_id}|{hit_id}|{hit_start}-{hit_end}")


def _group_n_terminal_hits(hits: list[Hit]) -> list[Hit]:
    """ Merge adjacent (<60 residues apart) AMP-binding N-terminal hits for one sequence into one hit each. """
    n_terminal = sorted((h for h in hits if h[0] == "AMP-binding"), key=lambda h: h[1])
    c_terminal = [h for h in hits if h[0] == "AMP-binding_C"]

    grouped: list[list[Hit]] = []
    if n_terminal:
        group = [n_terminal[0]]
        for i, hit_1 in enumerate(n_terminal):
            if i + 1 < len(n_terminal):
                hit_2 = n_terminal[i + 1]
                if hit_2[1] - hit_1[2] < 60:
                    group.append(hit_2)
                else:
                    grouped.append(group)
                    group = [hit_2]
            else:
                grouped.append(group)

    return [_merge_hits(g) for g in grouped] + c_terminal


def _hits_to_domains(
    id_to_hit: dict[str, HSP],
    path_in_fasta_file: str | Path,
    path_temp_dir: str | Path,
    use_profile_alignment: bool,
    hmm_version: int,
) -> list[AdenylationDomain]:
    hits_by_seq_id: dict[str, list[Hit]] = {}
    for hit_key in id_to_hit:
        seq_id, hit_id, hit_location = hit_key.split("|")
        hit_start_str, hit_end_str = hit_location.split("-")
        hits_by_seq_id.setdefault(seq_id, []).append((hit_id, int(hit_start_str), int(hit_end_str), hit_key))

    seq_id_to_domains: dict[str, list[AdenylationDomain]] = {}
    for seq_id, hits in hits_by_seq_id.items():
        grouped_hits = _group_n_terminal_hits(hits)

        for hit_id_1, hit_start_1, hit_end_1, hit_key_1 in grouped_hits:
            if hit_id_1 != "AMP-binding":
                continue

            seq_id_to_domains.setdefault(seq_id, [])
            match_found = False
            for hit_id_2, hit_start_2, hit_end_2, hit_key_2 in grouped_hits:
                if hit_id_2 == "AMP-binding_C" and hit_start_2 > hit_end_1 and (hit_start_2 - hit_end_1) < 200:
                    a_domain = AdenylationDomain(protein_name=seq_id, start=hit_start_1, end=hit_end_2)
                    if hmm_version == 2 and not use_profile_alignment:
                        a_domain.set_signatures_from_hmm(id_to_hit[hit_key_1], id_to_hit[hit_key_2])
                    seq_id_to_domains[seq_id].append(a_domain)
                    match_found = True
                    break

            if not match_found:
                a_domain = AdenylationDomain(protein_name=seq_id, start=hit_start_1, end=hit_end_1)
                if hmm_version == 2 and not use_profile_alignment:
                    a_domain.set_signatures_from_hmm(id_to_hit[hit_key_1])
                seq_id_to_domains[seq_id].append(a_domain)

    for domains in seq_id_to_domains.values():
        domains.sort(key=lambda d: d.start)

    fasta = parse_fasta_file(path_in_fasta_file)
    for seq_id, sequence in fasta.items():
        counter = 1
        for a_domain in seq_id_to_domains.get(seq_id, []):
            if seq_id != a_domain.protein_name:
                raise ValueError("protein name mismatch")
            domain_sequence = sequence[a_domain.start:a_domain.end]
            if len(domain_sequence) > 100:
                a_domain.sequence = domain_sequence
                a_domain.protein_sequence = sequence
                a_domain.domain_nr = counter
                counter += 1

    domains_out: list[AdenylationDomain] = []
    for a_domains in seq_id_to_domains.values():
        for a_domain in a_domains:
            if hmm_version == 2 and not use_profile_alignment:
                if a_domain.sequence and a_domain.extended_signature and a_domain.signature and a_domain.domain_nr:
                    domains_out.append(a_domain)
            elif use_profile_alignment:
                if not a_domain.sequence:
                    continue
                a_domain.set_signatures_from_profile(path_temp_dir)
                if a_domain.sequence and a_domain.extended_signature and a_domain.signature and a_domain.domain_nr:
                    domains_out.append(a_domain)
            else:
                if a_domain.sequence and a_domain.domain_nr:
                    domains_out.append(a_domain)

    domains_out.sort(key=lambda d: (d.protein_name, d.start))
    return domains_out


def _renumber_domains(domains: list[AdenylationDomain]) -> None:
    domains.sort(key=lambda d: (d.protein_name, d.start))
    protein_name = None
    counter = 0
    for domain in domains:
        counter = 1 if domain.protein_name != protein_name else counter + 1
        protein_name = domain.protein_name
        domain.domain_nr = counter


def _hmmer3_only_domains(
    hmmer2_domains: list[AdenylationDomain],
    hmmer3_domains: list[AdenylationDomain],
    path_temp_dir: str | Path,
) -> list[AdenylationDomain]:
    """ HMMER3 hits with no overlapping HMMER2 hit: extract their signatures via MUSCLE3 profile alignment. """
    unique: list[AdenylationDomain] = []
    for domain_3 in hmmer3_domains:
        if not any(
            domain_3.protein_name == domain_2.protein_name and domain_3.overlaps(domain_2, threshold=50)
            for domain_2 in hmmer2_domains
        ):
            domain_3.set_signatures_from_profile(path_temp_dir)
            unique.append(domain_3)
    return unique


def _extend_hmmer2_domains_with_hmmer3_bounds(
    hmmer2_domains: list[AdenylationDomain],
    hmmer3_domains: list[AdenylationDomain],
    path_in_fasta_file: str | Path,
) -> None:
    """ Where HMMER3 found a longer overlapping domain than HMMER2, widen the HMMER2 domain to match. """
    fasta = parse_fasta_file(path_in_fasta_file)
    for domain_2 in hmmer2_domains:
        for domain_3 in hmmer3_domains:
            if domain_2.protein_name == domain_3.protein_name and domain_2.overlaps(domain_3, threshold=50):
                if len(domain_3.sequence) > len(domain_2.sequence):
                    domain_2.start, domain_2.end = domain_3.start, domain_3.end
                    domain_2.sequence = fasta[domain_2.protein_name][domain_2.start:domain_2.end]


def get_domains(
    path_in_fasta_file: str | Path,
    path_temp_dir: str | Path,
    hmm_timeout: float | None = None,
    use_muscle_fallback: bool = True,
) -> list[AdenylationDomain]:
    """
    Extract every adenylation domain from a (protein) fasta file, following PARAS' default pipeline:
    HMMER2 signature extraction first, HMMER3 + MUSCLE3 profile alignment to pick up anything HMMER2 missed.

    :param path_in_fasta_file: protein fasta file to scan.
    :param path_temp_dir: scratch directory for intermediate HMMER/MUSCLE output.
    :param hmm_timeout: seconds to allow the hmmpfam2/hmmscan calls, or None for their own
        default -- override this for a large batch fasta file (e.g. retromol_paras.train's
        ~3654-domain training set), where the default (sized for a single/handful of domains)
        would time out a legitimately-slow full-file scan.
    :param use_muscle_fallback: whether to run a MUSCLE3 profile alignment (one subprocess
        call per affected domain, against a ~1000-sequence reference alignment -- the
        dominant per-domain cost in practice) to recover signatures for domains HMMER3
        detects but HMMER2's smaller/older profile misses entirely. False skips this
        step -- those domains are dropped (no signature, no prediction) instead of
        recovered, trading some coverage for speed.
    :return: AdenylationDomain instances, each with `signature`/`extended_signature` populated.
    """
    path_temp_dir = Path(path_temp_dir)
    path_temp_dir.mkdir(parents=True, exist_ok=True)
    ensure_hmm3_pressed(HMM3_FILE)

    hmm2_out = path_temp_dir / "run.hmm2_result"
    hmm3_out = path_temp_dir / "run.hmm3_result"
    hmm_kwargs = {} if hmm_timeout is None else {"timeout": hmm_timeout}
    run_hmmpfam2(HMM2_FILE, path_in_fasta_file, hmm2_out, **hmm_kwargs)
    run_hmmscan(HMM3_FILE, path_in_fasta_file, hmm3_out, **hmm_kwargs)

    id_to_hit_2 = parse_hmm_results(hmm2_out, hmmer_version=2)
    id_to_hit_3 = parse_hmm_results(hmm3_out, hmmer_version=3)

    domains_2 = _hits_to_domains(id_to_hit_2, path_in_fasta_file, path_temp_dir, use_profile_alignment=False, hmm_version=2)

    if use_muscle_fallback:
        domains_3 = _hits_to_domains(id_to_hit_3, path_in_fasta_file, path_temp_dir, use_profile_alignment=False, hmm_version=3)
        _extend_hmmer2_domains_with_hmmer3_bounds(domains_2, domains_3, path_in_fasta_file)
        unique_3 = _hmmer3_only_domains(domains_2, domains_3, path_temp_dir)
        domains = domains_2 + unique_3
    else:
        domains = domains_2

    _renumber_domains(domains)
    return domains


def extract_domains(
    path_in_fasta_file: str | Path,
    path_temp_dir: str | Path,
    hmm_timeout: float | None = None,
    use_muscle_fallback: bool = True,
) -> list[AdenylationDomain]:
    """
    Like `get_domains`, but safe for arbitrary fasta headers: renames sequences to plain
    integers before running HMMER (some header formats break HMMER's parser), then restores
    the original headers as `protein_name` afterwards. This is the entry point to use for a
    fasta file with headers you don't control.

    :param hmm_timeout: see `get_domains`.
    :param use_muscle_fallback: see `get_domains`.
    """
    path_temp_dir = Path(path_temp_dir)
    mapping_file, renamed_fasta_file = rename_sequences(path_in_fasta_file, path_temp_dir)
    domains = get_domains(renamed_fasta_file, path_temp_dir, hmm_timeout=hmm_timeout, use_muscle_fallback=use_muscle_fallback)

    new_to_original = load_renaming_map(mapping_file)
    for domain in domains:
        domain.protein_name = new_to_original[domain.protein_name]

    return domains
