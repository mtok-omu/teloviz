"""Tests for the HTML telomere report (spec section 10)."""

from __future__ import annotations

from teloviz.calling import Call
from teloviz.report import write_report

_CALLS = [
    Call("chr1", 100_000, 300, 250, True, True),    # both
    Call("chr2", 80_000, 300, 0, True, False),      # 5' only
    Call("chr3", 60_000, 0, 0, False, False),       # none
]
_META = {
    "command": "teloviz in.tsv --min-count 50 -o out",
    "input": "in.tsv", "mode": "both", "dist_kb": "30", "call_min": 50,
    "min_count": 50, "motif": None, "window_size": 10000,
}


def test_report_written_with_summary_and_settings(tmp_path):
    p = write_report(tmp_path / "r.html", _CALLS, meta=_META)
    assert p.exists()
    html = p.read_text()
    # Summary counts.
    assert "both ends capped: <b>1</b>" in html
    assert "5' capped: 2" in html and "3' capped: 1" in html
    # Settings echo the command and thresholds.
    assert "teloviz in.tsv --min-count 50 -o out" in html
    assert "30 kb" in html and "50" in html


def test_report_rows_and_status(tmp_path):
    html = write_report(tmp_path / "r.html", _CALLS, meta=_META).read_text()
    assert "chr1 *" in html          # both-ends chromosome gets an asterisk
    assert "chr2" in html and "chr3" in html
    assert "5&#x27; only" in html    # 5'-only status (apostrophe escaped)
    assert ">none<" in html


def test_report_creates_parent_dirs(tmp_path):
    out = tmp_path / "nested" / "sub" / "r.html"
    p = write_report(out, _CALLS, meta=_META)
    assert p.exists()
