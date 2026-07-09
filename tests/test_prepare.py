"""Tests for the data-preparation stage (spec section 4)."""

from __future__ import annotations

import pandas as pd

from teloviz.prepare import DEFAULT_WINDOW_SIZE, estimate_window_size, prepare


def _df(rows):
    return pd.DataFrame(rows, columns=["id", "window", "forward", "reverse", "motif"])


def test_window_size_estimated_once_over_whole_file():
    # chr1 spaced by 10000; chr2 has a single window (no diff of its own).
    df = _df([
        ("chr1", 10000, 1, 0, "T"),
        ("chr1", 20000, 1, 0, "T"),
        ("chr1", 30000, 1, 0, "T"),
        ("chr2", 10000, 1, 0, "T"),
    ])
    assert estimate_window_size(df) == 10000


def test_window_size_fallback_when_all_single_window():
    df = _df([("chr1", 10000, 1, 0, "T"), ("chr2", 50000, 1, 0, "T")])
    assert estimate_window_size(df) == DEFAULT_WINDOW_SIZE


def test_motifs_aggregated_per_window():
    df = _df([
        ("chr1", 10000, 2, 3, "TTAGGG"),
        ("chr1", 10000, 5, 1, "TTTAGGG"),
    ])
    p = prepare(df)
    row = p.table[(p.table.id == "chr1") & (p.table.window == 10000)].iloc[0]
    assert row["forward"] == 7 and row["reverse"] == 4


def test_natural_sort_and_min_len(tmp_path):
    df = _df([
        ("chr1", 10000, 1, 0, "T"),
        ("chr2", 10000, 1, 0, "T"),
        ("chr10", 10000, 1, 0, "T"),
        ("scaf", 5000, 1, 0, "T"),
    ])
    # .fai gives chr* = 100000, scaf = 5000; drop scaf with min_len.
    fai = {"chr1": 100000, "chr2": 100000, "chr10": 100000, "scaf": 5000}
    p = prepare(df, fai=fai, min_len=10000)
    assert p.order == ["chr1", "chr2", "chr10"]  # chr2 before chr10 (natural)
    assert "scaf" not in p.lengths


def test_length_from_fai_else_max_window():
    df = _df([("chr1", 10000, 1, 0, "T"), ("chr1", 20000, 1, 0, "T")])
    assert prepare(df).lengths["chr1"] == 20000
    assert prepare(df, fai={"chr1": 33333}).lengths["chr1"] == 33333
