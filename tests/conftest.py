"""Shared fixtures: in-process CLI runner and temporary input files."""
import io
import sys

import pytest

from fastxkit import cli


@pytest.fixture
def run_cli(monkeypatch, capsys):
    """Returns run(*argv, stdin=b'') -> (exit_code, stdout, stderr)."""

    def run(*argv, stdin=b""):
        reader = io.BufferedReader(io.BytesIO(stdin))
        monkeypatch.setattr(sys, "stdin", io.TextIOWrapper(reader))
        code = cli.main(list(argv))
        captured = capsys.readouterr()
        return code, captured.out, captured.err

    return run


@pytest.fixture
def make_file(tmp_path):
    """Returns make(name, content) -> path string; content is str or bytes."""

    def make(name, content):
        path = tmp_path / name
        if isinstance(content, str):
            path.write_text(content, newline="")
        else:
            path.write_bytes(content)
        return str(path)

    return make
