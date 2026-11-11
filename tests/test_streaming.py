"""Regression tests for actionable batch parsing errors."""

from retromol.io.streaming import _process_compound


def test_worker_error_includes_input_smiles_and_exception_type() -> None:
    result, error = _process_compound(("not-a-smiles", {"id": "broken row"}))

    assert result is None
    assert error is not None
    assert "ValueError" in error
    assert "SMILES='not-a-smiles'" in error
