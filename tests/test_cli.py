"""Smoke tests for the teloviz CLI skeleton."""

from __future__ import annotations

import pytest

from teloviz import __version__
from teloviz.cli import build_parser, main, parse_formats

WIN_TSV = "id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\nchr1\t10000\t2\t344\tTTAGGG\n"


def test_version_flag(capsys):
    with pytest.raises(SystemExit) as exc:
        build_parser().parse_args(["--version"])
    assert exc.value.code == 0
    assert __version__ in capsys.readouterr().out


def test_missing_input_returns_2(capsys):
    code = main(["does_not_exist.tsv"])
    assert code == 2
    assert "input not found" in capsys.readouterr().err


def test_parse_formats_rejects_unknown():
    with pytest.raises(Exception):
        parse_formats("pdf,jpeg")


def test_parse_formats_accepts_pdf_png_svg():
    assert parse_formats("pdf, png ,svg") == ["pdf", "png", "svg"]


def test_out_prefix_derived_from_input(tmp_path, capsys):
    tsv = tmp_path / "sample_telomeric_repeat_windows.tsv"
    tsv.write_text(WIN_TSV)
    code = main([str(tsv), "--mode", "both"])
    assert code == 0
    out = capsys.readouterr().out
    assert "out_prefix=sample" in out
    assert "mode=both" in out
    assert "formats=pdf" in out
    assert "parsed OK" in out
