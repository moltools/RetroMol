"""Lossless, streaming result storage in plain JSON and gzip JSONL."""

import gzip
import json
from pathlib import Path

import pytest

from retromol.io.json import dumps_compact, iter_json, open_json_file
from retromol.io.streaming import ResultEvent, stream_json_records
from retromol.model.result import Result
from retromol.model.submission import Submission
from retromol.pipelines.parsing import run_retromol
from retromol_cli import cli as cli_module


RECORDS = [
    {'id': 'same', 'note': 'Ångström, α/β; keep spaces inside strings', 'nested': {'flag': False, 'missing': None}},
    {'id': 'same', 'note': 'line one\nline two', 'numbers': [0, -3, 1.25]},
]


@pytest.mark.parametrize('suffix', ['.jsonl', '.jsonl.gz', '.ndjson', '.ndjson.gz', '.json', '.json.gz'])
def test_streaming_plain_and_gzip_preserve_all_json_values(tmp_path, suffix):
    path = tmp_path / ('records' + suffix)
    with open_json_file(path, 'wt') as handle:
        if suffix.startswith('.jsonl') or suffix.startswith('.ndjson'):
            for record in RECORDS:
                handle.write(dumps_compact(record) + '\n')
        else:
            handle.write(dumps_compact(RECORDS))
    assert list(iter_json(path)) == RECORDS
    assert list(stream_json_records(str(path))) == RECORDS
    if suffix.endswith('.gz'):
        assert path.read_bytes().startswith(b'\x1f\x8b')


def test_compact_json_only_removes_formatting_whitespace():
    text = dumps_compact(RECORDS)
    assert json.loads(text) == RECORDS
    assert 'Ångström' in text
    assert len(text.encode()) < len(json.dumps(RECORDS).encode())


def test_appending_gzip_runs_reads_all_members(tmp_path):
    for record in RECORDS:
        handle, path = cli_module._open_jsonl(str(tmp_path), None)
        with handle:
            handle.write(dumps_compact(record) + '\n')
    assert path == str(tmp_path / 'results.jsonl.gz')
    assert list(iter_json(path)) == RECORDS


def test_custom_jsonl_filename_supports_current_directory(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    handle, path = cli_module._open_jsonl(str(tmp_path), 'custom.jsonl')
    with handle:
        handle.write(dumps_compact(RECORDS[0]) + '\n')
    assert path == 'custom.jsonl.gz'
    assert list(iter_json(path)) == RECORDS[:1]


@pytest.fixture
def batch_run(monkeypatch, tmp_path, ruleset):
    result = run_retromol(Submission('CCO.[Na+]', props=RECORDS[0]), ruleset)
    payload = result.to_dict()
    input_path = tmp_path / 'input.jsonl.gz'
    with open_json_file(input_path, 'wt') as handle:
        handle.write(dumps_compact({'smiles': 'CCO.[Na+]', **RECORDS[0]}) + '\n')
    output = tmp_path / 'output'
    monkeypatch.setattr(cli_module, 'setup_logging', lambda **kwargs: None)
    monkeypatch.setattr(cli_module, 'add_file_handler', lambda *args, **kwargs: None)

    def events(**kwargs):
        rows = list(kwargs['row_iter'])
        assert rows == [{'smiles': 'CCO.[Na+]', **RECORDS[0]}]
        yield ResultEvent(payload, None)
        yield ResultEvent(None, 'intentional failed input')
        yield ResultEvent(payload, None)  # repeated IDs must not overwrite results

    monkeypatch.setattr(cli_module, 'run_retromol_stream', events)

    def run(*options):
        monkeypatch.setattr('sys.argv', ['retromol', '-o', str(output), 'batch', '--json', str(input_path), '--no-tqdm', *options])
        cli_module.main()
        return output, payload

    return run


@pytest.mark.parametrize('compression,suffix', [('gzip', '.jsonl.gz'), ('none', '.jsonl')])
def test_batch_writes_lossless_results_with_selected_compression(batch_run, compression, suffix):
    output, payload = batch_run('--compression', compression)
    path = output / ('results' + suffix)
    records = list(iter_json(path))
    assert records == [payload, payload]
    for record in records:
        loaded = Result.from_dict(record)
        assert loaded.props == RECORDS[0]
        assert loaded.calculate_coverage() == .75
        assert len(loaded.reaction_rules_hash) == len(loaded.matching_rules_hash) == 64
    with open_json_file(path) as handle:
        assert handle.readline().rstrip('\n') == dumps_compact(payload)


def test_batch_compresses_by_default(batch_run):
    output, payload = batch_run()
    assert not (output / 'results.jsonl').exists()
    assert list(iter_json(output / 'results.jsonl.gz')) == [payload, payload]


@pytest.mark.parametrize('error', [RuntimeError, KeyboardInterrupt])
def test_interrupted_batch_closes_gzip_and_preserves_completed_records(batch_run, monkeypatch, error):
    original = cli_module.run_retromol_stream

    def interrupted(**kwargs):
        yield next(original(**kwargs))
        raise error('interrupted run')

    monkeypatch.setattr(cli_module, 'run_retromol_stream', interrupted)
    with pytest.raises(error):
        batch_run()
    path = Path(cli_module.cli().outdir) / 'results.jsonl.gz'
    records = list(iter_json(path))
    assert len(records) == 1
    assert Result.from_dict(records[0]).props == RECORDS[0]


@pytest.mark.parametrize('compression,suffix', [('gzip', '.json.gz'), ('none', '.json')])
def test_separate_result_files_preserve_repeated_ids(batch_run, compression, suffix):
    output, payload = batch_run('--results', 'files', '--compression', compression)
    paths = sorted(output.glob('result_*'))
    assert [p.name for p in paths] == [f'result_000000001{suffix}', f'result_000000003{suffix}']
    before = [p.read_bytes() for p in paths]
    for path in paths:
        with open_json_file(path) as handle:
            assert json.load(handle) == payload
    with pytest.raises(FileExistsError):
        batch_run('--results', 'files', '--compression', compression)
    assert [p.read_bytes() for p in paths] == before


def test_plain_output_rejects_misleading_gzip_extension(tmp_path, monkeypatch, capsys):
    output = tmp_path / 'output'
    monkeypatch.setattr('sys.argv', ['retromol', '-o', str(output), 'batch', '--json', 'input.jsonl',
                                    '--compression', 'none', '--jsonl-path', 'results.jsonl.gz'])
    with pytest.raises(SystemExit) as error:
        cli_module.main()
    assert error.value.code == 2
    assert 'requires --compression gzip' in capsys.readouterr().err
    assert not output.exists()


def test_single_can_save_compressed_result(tmp_path, monkeypatch):
    monkeypatch.setattr(cli_module, 'setup_logging', lambda **kwargs: None)
    monkeypatch.setattr(cli_module, 'add_file_handler', lambda *args, **kwargs: None)
    monkeypatch.setattr(cli_module, 'visualize_reaction_graph', lambda *args, **kwargs: None)
    monkeypatch.setattr('sys.argv', ['retromol', '-o', str(tmp_path), 'single', '-s', 'CCO.[Na+]',
                                    '--compression', 'gzip', '--max-assembly-graphs', '0'])
    cli_module.main()
    assert not (tmp_path / 'result.json').exists()
    with gzip.open(tmp_path / 'result.json.gz', 'rt') as handle:
        result = Result.from_dict(json.load(handle))
    assert result.calculate_coverage() == .75
