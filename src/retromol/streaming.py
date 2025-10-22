"""Streaming RetroMol runs with multiprocessing."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from multiprocessing import Pool
from typing import Any

import yaml
from pandas import DataFrame, read_csv
from rdkit.Chem.rdmolfiles import SDMolSupplier

from retromol import api, io, rules
from retromol.chem import mol_to_inchikey, mol_to_smiles, sanitize_mol
from retromol.helpers import iter_json
from retromol.rules import Rules

_G_RULE_SET = None
_G_WAVE_CONFIGS = None
_G_MATCH_STEREO = None


def _init_worker(rule_set: Rules, wave_configs: list[dict[str, Any]] | None, match_stereo: bool) -> None:
    """
    Initialize worker process with necessary global variables.

    :param rule_set: reaction/matching rule set
    :param wave_configs: wave configuration dicts
    :param match_stereo: whether to match stereo in RetroMol runs
    """
    global _G_RULE_SET, _G_WAVE_CONFIGS, _G_MATCH_STEREO
    _G_RULE_SET = rule_set
    _G_WAVE_CONFIGS = wave_configs
    _G_MATCH_STEREO = match_stereo


def _process_compound(
    args_tuple: tuple[str, str, dict[str, Any]],
) -> tuple[str, dict[str, Any] | None, str | None]:
    """
    Process a single compound in a worker process.

    :param args_tuple: (inchikey, smiles, props)
    :return: (inchikey, serialized_result or None on error, error message or None on success)
    """
    inchikey, smiles, props = args_tuple
    try:
        mol = io.Input(inchikey, smiles, props=props or {})
        if _G_RULE_SET is None:
            raise RuntimeError("Worker not properly initialized with rule set.")
        if _G_WAVE_CONFIGS is None:
            raise RuntimeError("Worker not properly initialized with wave configs.")
        result_obj = api.run_retromol_with_timeout(
            mol,
            _G_RULE_SET,
            _G_WAVE_CONFIGS,
            _G_MATCH_STEREO if _G_MATCH_STEREO is not None else False,
        )
        return inchikey, result_obj.serialize(), None
    except Exception as e:
        # traceback not returned here to keep workers light-weight; caller can log
        return inchikey, None, str(e)


@dataclass
class ResultEvent:
    """
    Represents the result of processing a single compound.

    :param inchikey: InChIKey of the processed compound
    :param result: serialized result dict or None if there was an error
    :param error: error message string or None if processing was successful
    """

    inchikey: str
    result: dict[str, Any] | None  # serialized result or None on error
    error: str | None  # error message or None on success


def _task_buffered_iterator(
    source_iter: Iterable[dict[str, Any]],
    *,
    id_col: str,
    smiles_col: str,
    batch_size: int,
) -> Iterator[list[tuple[str, str, dict[str, Any]]]]:
    """
    Convert row dicts into (inchikey, smiles, props) tuples and yield in batches.

    :param source_iter: iterable of row dicts
    :param id_col: name of column containing InChIKey
    :param smiles_col: name of column containing SMILES
    :param batch_size: number of compounds per batch
    :return: iterator over lists of (inchikey, smiles, props) tuples
    """
    buf: list[tuple[str, str, dict[str, Any]]] = []
    for rec in source_iter:
        if id_col not in rec or smiles_col not in rec:
            continue
        ik = str(rec[id_col])
        smi = str(rec[smiles_col])
        buf.append((ik, smi, rec))
        if len(buf) >= batch_size:
            yield buf
            buf = []
    if buf:
        yield buf


def run_retromol_stream(
    *,
    # Either provide loaded objects...
    rule_set: Rules | None = None,
    wave_configs: list[dict[str, Any]] | None = None,
    # ...or point to files (if provided, they override the loaded objects)
    reaction_rules_path: str | None = None,
    matching_rules_path: str | None = None,
    wave_config_path: str | None = None,
    match_stereo: bool = False,
    # Data source: an iterable of row dicts containing id_col and smiles_col
    row_iter: Iterable[dict[str, Any]],
    id_col: str = "inchikey",
    smiles_col: str = "smiles",
    # Concurrency knobs (match CLI defaults)
    workers: int = 1,
    batch_size: int = 2000,
    pool_chunksize: int = 50,
    maxtasksperchild: int = 2000,
    # Optional sink callback
    on_result: Callable[[ResultEvent], None] | None = None,
) -> Iterator[ResultEvent]:
    """
    Stream RetroMol results with multiprocessing, yielding ResultEvent as soon as
    each compound finishes. No files/logs are written here—callers are free to do so.

    :param rule_set: pre-loaded reaction/matching rule set
    :param wave_configs: pre-loaded wave configuration dicts
    :param reaction_rules_path: path to reaction rules file (YAML)
    :param matching_rules_path: path to matching rules file (YAML)
    :param wave_config_path: path to wave configuration file (YAML)
    :param match_stereo: whether to match stereo in RetroMol runs
    :param row_iter: iterable of row dicts containing at least id_col and smiles_col
    :param id_col: name of column containing InChIKey (default: "inchikey")
    :param smiles_col: name of column containing SMILES (default: "smiles")
    :param workers: number of worker processes (default: 1)
    :param batch_size: number of compounds to send to each worker at once (default: 2000)
    :param pool_chunksize: chunksize for imap_unordered (default: 50)
    :param maxtasksperchild: max tasks per worker before restart (default: 2000)
    :param on_result: optional callback receiving each ResultEvent as it arrives
    :return: iterator over ResultEvent objects
    """
    # Load/prepare config exactly once (like CLI)
    if reaction_rules_path and matching_rules_path:
        rule_set = rules.load_rules_from_files(reaction_rules_path, matching_rules_path)
    elif rule_set is None:
        raise ValueError("Provide either (reaction_rules_path & matching_rules_path) or an already loaded rule_set.")

    if wave_config_path:
        with open(wave_config_path) as f:
            wave_configs = yaml.safe_load(f)
    if wave_configs is None:
        raise ValueError("Provide wave_configs dict or wave_config_path.")

    # Start worker pool with same init pattern
    with Pool(
        processes=workers,
        initializer=_init_worker,
        initargs=(rule_set, wave_configs, match_stereo),
        maxtasksperchild=maxtasksperchild,
    ) as pool:
        for task_batch in _task_buffered_iterator(
            row_iter, id_col=id_col, smiles_col=smiles_col, batch_size=batch_size
        ):
            for ik, serialized, err in pool.imap_unordered(_process_compound, task_batch, chunksize=pool_chunksize):
                evt = ResultEvent(ik, serialized, err)
                if on_result is not None:
                    on_result(evt)
                yield evt


def stream_table_rows(
    path: str,
    *,
    sep: str = ",",
    chunksize: int = 20_000,
) -> Iterator[dict[str, Any]]:
    """
    Stream CSV/TSV rows as dicts. Keeps memory usage low (chunked).

    :param path: path to CSV/TSV file
    :param sep: field separator (default: ",")
    :param chunksize: number of rows to read per chunk (default: 20,000)
    :return: iterator over row dicts
    """
    chunks: Iterator[DataFrame] = read_csv(
        path,
        sep=sep,
        chunksize=chunksize,
        dtype=str,
        keep_default_na=False,
    )

    for chunk in chunks:
        # iterrows() -> Iterator[Tuple[int, Series]]
        for _, row in chunk.iterrows():
            yield row.to_dict()


def stream_sdf_records(
    sdf_path: str,
    *,
    fast: bool = False,
) -> Iterator[dict[str, Any]]:
    """
    Stream SDF as dict rows: {'inchikey': <IK>, 'smiles': <SMI>, ...props}
    Matches CLI behavior including opportunistic sanitize for IK/SMILES.

    :param sdf_path: path to SDF file
    :param fast: if True, skips sanitization and H removal (default: False)
    :return: iterator over record dicts
    """
    sanitize = not fast
    removeHs = fast
    suppl = SDMolSupplier(sdf_path, sanitize=sanitize, removeHs=removeHs)
    for mol in suppl:
        if mol is None:
            continue
        try:
            try:
                ik = mol_to_inchikey(mol)
                smi = mol_to_smiles(mol)
            except Exception:
                sanitize_mol(mol)
                ik = mol_to_inchikey(mol)
                smi = mol_to_smiles(mol)
            rec = {"inchikey": ik, "smiles": smi}
            for pname in mol.GetPropNames():
                rec[pname] = mol.GetProp(pname)
            yield rec
        except Exception:
            continue


def stream_json_records(
    path: str,
    *,
    jsonl: bool = False,
) -> Iterator[dict[str, Any]]:
    """
    Stream JSON or JSONL records as dicts.

    :param path: path to JSON or JSONL file
    :param jsonl: if True, treat as JSONL (one JSON object per line)
    :return: iterator over record dicts
    """
    for rec in iter_json(path, jsonl=jsonl):
        if isinstance(rec, dict):
            yield rec
