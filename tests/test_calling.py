"""Tests for telomere-cap calling (spec section 10)."""

from __future__ import annotations

import pandas as pd

from teloviz.calling import call_telomeres
from teloviz.prepare import prepare


def _prepared(rows, fai=None):
    df = pd.DataFrame(rows, columns=["id", "window", "forward", "reverse", "motif"])
    return prepare(df, fai=fai)


def test_end_regions_summed_independently():
    # chr1 length 100kb; 5' has 60 within 30kb, 3' has 10 within 30kb.
    rows = [("chr1", w, 0, 0, "T") for w in range(10000, 100001, 10000)]
    rows[0] = ("chr1", 10000, 40, 20, "T")     # 5' window: total 60
    rows[-1] = ("chr1", 100000, 6, 4, "T")     # 3' window: total 10
    p = _prepared(rows, fai={"chr1": 100000})
    (c,) = call_telomeres(p, dist_bp=30000, threshold=50)
    assert c.five_count == 60 and c.three_count == 10
    assert c.five and not c.three
    assert c.status == "5' only" and not c.both


def test_threshold_is_inclusive():
    rows = [("chr1", 10000, 50, 0, "T"), ("chr1", 20000, 0, 0, "T")]
    (c,) = call_telomeres(_prepared(rows, fai={"chr1": 20000}), dist_bp=10000, threshold=50)
    assert c.five  # exactly at threshold counts as capped


def test_dist_controls_how_far_in_we_look():
    # A telomere sitting at 20-30kb is seen at dist=30kb but not at dist=10kb.
    rows = [("chr1", 10000, 0, 0, "T"), ("chr1", 30000, 200, 0, "T"),
            ("chr1", 60000, 0, 0, "T")]
    p = _prepared(rows, fai={"chr1": 60000})
    assert not call_telomeres(p, dist_bp=10000, threshold=50)[0].five
    assert call_telomeres(p, dist_bp=30000, threshold=50)[0].five


def test_both_ends_capped_status():
    rows = [("chr1", 10000, 300, 0, "T"), ("chr1", 50000, 0, 0, "T"),
            ("chr1", 90000, 0, 300, "T")]
    (c,) = call_telomeres(_prepared(rows, fai={"chr1": 90000}), dist_bp=30000, threshold=50)
    assert c.both and c.status == "both"


def test_uses_min_count_filtered_counts():
    # Windows below --min-count are zeroed in prepare(), so they don't count
    # toward the call either (calling shares the plotted table).
    rows = [("chr1", 10000, 30, 10, "T"),   # total 40 -> zeroed by min_count=50
            ("chr1", 20000, 0, 0, "T")]
    df = pd.DataFrame(rows, columns=["id", "window", "forward", "reverse", "motif"])
    p = prepare(df, fai={"chr1": 20000}, min_count=50)
    (c,) = call_telomeres(p, dist_bp=20000, threshold=10)
    assert c.five_count == 0 and not c.five
