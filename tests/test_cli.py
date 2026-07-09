"""Smoke + end-to-end tests for the teloviz CLI."""

from __future__ import annotations

import pytest

from teloviz import __version__
from teloviz.cli import build_parser, main, parse_formats

HEADER = "id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\n"


def _windows_tsv() -> str:
    rows = []
    for cid, length in (("chr1", 5), ("chr2", 4), ("chr10", 6)):
        for k in range(1, length + 1):
            w = k * 10000
            fwd = 300 if k == 1 else (0 if k < length else 2)
            rev = 2 if k == 1 else (0 if k < length else 300)
            rows.append(f"{cid}\t{w}\t{fwd}\t{rev}\tTTAGGG")
    return HEADER + "\n".join(rows) + "\n"


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


def test_end_to_end_both_modes(tmp_path):
    tsv = tmp_path / "sample_telomeric_repeat_windows.tsv"
    tsv.write_text(_windows_tsv())
    prefix = tmp_path / "out"
    code = main([str(tsv), "--mode", "both", "--format", "pdf,png", "-o", str(prefix)])
    assert code == 0
    for mode in ("sum", "orientation"):
        for fmt in ("pdf", "png"):
            f = tmp_path / f"out.{mode}.{fmt}"
            assert f.exists() and f.stat().st_size > 0


def test_min_len_drops_all(tmp_path, capsys):
    tsv = tmp_path / "x_telomeric_repeat_windows.tsv"
    tsv.write_text(_windows_tsv())
    code = main([str(tsv), "--min-len", "10_000_000", "-o", str(tmp_path / "out")])
    assert code == 1
    assert "no sequences" in capsys.readouterr().err
