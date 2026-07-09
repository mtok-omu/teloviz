"""Smoke tests for the teloviz CLI skeleton."""

from __future__ import annotations

import pytest

from teloviz import __version__
from teloviz.cli import DEFAULT_MOTIF, build_parser, main


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        build_parser().parse_args(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_missing_input_returns_2(capsys):
    code = main(["--input", "does_not_exist.fasta"])
    assert code == 2
    assert "input not found" in capsys.readouterr().err


def test_parse_ok(tmp_path, capsys):
    fa = tmp_path / "g.fasta"
    fa.write_text(">chr1\nACGT\n")
    code = main(["--input", str(fa)])
    assert code == 0
    out = capsys.readouterr().out
    assert DEFAULT_MOTIF in out
    assert "parsed OK" in out
