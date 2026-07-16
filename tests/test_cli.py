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


def test_proximity_default_is_500kb():
    # Spec section 5.4: a real NOR can sit well inside an end (often >100 kb).
    args = build_parser().parse_args(["x.tsv"])
    assert args.proximity == 500.0


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


def test_writes_telomere_report_by_default(tmp_path):
    tsv = tmp_path / "sample_telomeric_repeat_windows.tsv"
    tsv.write_text(_windows_tsv())
    prefix = tmp_path / "out"
    code = main([str(tsv), "-o", str(prefix), "--format", "png"])
    assert code == 0
    report = tmp_path / "out.telomere_report.html"
    assert report.exists() and report.stat().st_size > 0
    assert "telomere-cap report" in report.read_text()


def test_no_call_suppresses_report(tmp_path):
    tsv = tmp_path / "x_telomeric_repeat_windows.tsv"
    tsv.write_text(_windows_tsv())
    prefix = tmp_path / "out"
    code = main([str(tsv), "-o", str(prefix), "--format", "png", "--no-call"])
    assert code == 0
    assert not (tmp_path / "out.telomere_report.html").exists()


def _features_bed() -> str:
    # chr10 (length 60000) 3' end: a 45S array right at the end.
    return "chr10\t54000\t59000\t45S_array\trdna_45S\n"


def test_rdna_adds_feature_notes_to_report(tmp_path):
    tsv = tmp_path / "sample_telomeric_repeat_windows.tsv"
    tsv.write_text(_windows_tsv())
    bed = tmp_path / "features.bed"
    bed.write_text(_features_bed())
    prefix = tmp_path / "out"
    # Force every end un-capped so the 3'-end feature is actually annotated.
    code = main([str(tsv), "--rDNA", str(bed), "--call-min", "100000",
                 "-o", str(prefix), "--format", "png"])
    assert code == 0
    html = (tmp_path / "out.telomere_report.html").read_text()
    assert "feature notes" in html          # feature column present
    assert "rdna_45S (45S_array) 1 kb from the 3' end" in html


def test_default_has_no_feature_column(tmp_path):
    tsv = tmp_path / "sample_telomeric_repeat_windows.tsv"
    tsv.write_text(_windows_tsv())
    prefix = tmp_path / "out"
    code = main([str(tsv), "-o", str(prefix), "--format", "png"])
    assert code == 0
    assert "feature notes" not in (tmp_path / "out.telomere_report.html").read_text()


def test_missing_rdna_bed_returns_2(tmp_path, capsys):
    tsv = tmp_path / "sample_telomeric_repeat_windows.tsv"
    tsv.write_text(_windows_tsv())
    code = main([str(tsv), "--rDNA", "no_such.bed", "-o", str(tmp_path / "out")])
    assert code == 2
    assert "feature BED not found" in capsys.readouterr().err


def test_min_len_drops_all(tmp_path, capsys):
    tsv = tmp_path / "x_telomeric_repeat_windows.tsv"
    tsv.write_text(_windows_tsv())
    code = main([str(tsv), "--min-len", "10_000_000", "-o", str(tmp_path / "out")])
    assert code == 1
    assert "no sequences" in capsys.readouterr().err
