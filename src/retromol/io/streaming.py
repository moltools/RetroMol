"""Streaming RetroMol runs with multiprocessing."""

from __future__ import annotations

from dataclasses import dataclass
from multiprocessing import Pool
from typing import Any, Callable, Iterable, Iterator

from pandas import DataFrame, read_csv
from rdkit.Chem.rdmolfiles import SDMolSupplier

from retromol.model.rules import RuleSet
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol_with_timeout
from retromol.chem.mol import mol_to_smiles, sanitize_mol
from retromol.io.json import iter_json

_G_RULESET = None


def _init_worker(ruleset: RuleSet) -> None:
    """
    Initialize worker process with necessary global variables.

    :param ruleset: reaction/matching rule set
    :param wave_configs: wave configuration dicts
    """
    global _G_RULESET
    _G_RULESET = ruleset


def _process_compound(args_tuple: tuple[str, dict[str, Any]]) -> tuple[dict[str, Any] | None, str | None]:
    """
    Process a single compound in a worker process.

    :param args_tuple: (smiles, props)
    :return: (serialized_result or None on error, error message or None on success)
    """
    smiles, props = args_tuple
    try:
        submission = Submission(smiles, props=props)
        if _G_RULESET is None:
            raise RuntimeError("worker not properly initialized with rule set")
        result_obj = run_retromol_with_timeout(submission, _G_RULESET)
        return result_obj.to_dict(), None
    except Exception as e:
        # Traceback not returned here to keep workers light-weight; caller can log
        return None, str(e)


@dataclass
class ResultEvent:
    """
    Represents the result of processing a single compound.

    :param result: serialized result dict or None if there was an error
    :param error: error message string or None if processing was successful
    """

    result: dict[str, Any] | None  # serialized result or None on error
    error: str | None  # error message or None on success


def _task_buffered_iterator(
    source_iter: Iterable[dict[str, Any]],
    smiles_col: str,
    batch_size: int,
) -> Iterator[list[tuple[str, dict[str, Any]]]]:
    """
    Convert row dicts into (smiles, props) tuples and yield in batches.

    :param source_iter: iterable of row dicts
    :param smiles_col: name of column containing SMILES
    :param batch_size: number of compounds per batch
    :return: iterator over lists of (smiles, props) tuples
    """
    buf: list[tuple[str, dict[str, Any]]] = []
    for rec in source_iter:
        if smiles_col not in rec:
            continue
        smi = str(rec[smiles_col])
        buf.append((smi, rec))
        if len(buf) >= batch_size:
            yield buf
            buf = []
    if buf:
        yield buf


def run_retromol_stream(
    ruleset: RuleSet,
    row_iter: Iterable[dict[str, Any]],
    smiles_col: str = "smiles",
    workers: int = 1,
    batch_size: int = 2000,
    pool_chunksize: int = 50,
    maxtasksperchild: int = 2000,
    on_result: Callable[[ResultEvent], None] | None = None,
) -> Iterator[ResultEvent]:
    """
    Stream RetroMol results with multiprocessing, yielding ResultEvent as soon as
    each compound finishes. No files/logs are written here—callers are free to do so.

    :param ruleset: pre-loaded reaction/matching rule set
    :param row_iter: iterable of row dicts containing at least id_col and smiles_col
    :param smiles_col: name of column containing SMILES (default: "smiles")
    :param workers: number of worker processes (default: 1)
    :param batch_size: number of compounds to send to each worker at once (default: 2000)
    :param pool_chunksize: chunksize for imap_unordered (default: 50)
    :param maxtasksperchild: max tasks per worker before restart (default: 2000)
    :param on_result: optional callback receiving each ResultEvent as it arrives
    :return: iterator over ResultEvent objects
    """
    # Start worker pool with same init pattern
    with Pool(
        processes=workers,
        initializer=_init_worker,
        initargs=(ruleset,),
        maxtasksperchild=maxtasksperchild,
    ) as pool:
        for task_batch in _task_buffered_iterator(row_iter, smiles_col=smiles_col, batch_size=batch_size):
            for serialized, err in pool.imap_unordered(_process_compound, task_batch, chunksize=pool_chunksize):
                evt = ResultEvent(serialized, err)
                if on_result is not None:
                    on_result(evt)
                yield evt


def stream_table_rows(path: str, sep: str = ",", chunksize: int = 20_000) -> Iterator[dict[str, Any]]:
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


def stream_sdf_records(sdf_path: str, fast: bool = False) -> Iterator[dict[str, Any]]:
    """
    Stream SDF as dict rows: {'smiles': <SMI>, ...props}.

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
                smi = mol_to_smiles(mol)
            except Exception:
                sanitize_mol(mol, fix_hydrogens=False)
                smi = mol_to_smiles(mol)
            rec = {"smiles": smi}
            for pname in mol.GetPropNames():
                rec[pname] = mol.GetProp(pname)
            yield rec
        except Exception:
            continue


def stream_json_records(path: str, jsonl: bool = False) -> Iterator[dict[str, Any]]:
    """
    Stream JSON or JSONL records as dicts.

    :param path: path to JSON or JSONL file
    :param jsonl: if True, treat as JSONL (one JSON object per line)
    :return: iterator over record dicts
    """
    for rec in iter_json(path, jsonl=jsonl):
        if isinstance(rec, dict):
            yield rec
