"""Tests for the HTML telomere report (spec section 10)."""

from __future__ import annotations

from teloviz.calling import Call
from teloviz.features import Feature, FeatureSet
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


def test_report_has_no_feature_column_without_features(tmp_path):
    html = write_report(tmp_path / "r.html", _CALLS, meta=_META).read_text()
    assert "feature notes" not in html


# chr7 3' un-capped with a 45S array right at that end; chr8 both ends bare, no
# feature anywhere -> "no feature nearby".
_FEAT_CALLS = [
    Call("chr7", 700_000, 300, 0, True, False),   # 3' un-capped
    Call("chr8", 500_000, 0, 0, False, False),    # both un-capped, no feature
]
_FS = FeatureSet(
    features=[Feature("chr7", 690_000, 699_000, "45S_array_n28", "rdna_45S")],
    lane_order=["rdna_45S"],
)


def test_report_annotates_uncapped_end_with_nearby_feature(tmp_path):
    html = write_report(tmp_path / "r.html", _FEAT_CALLS, meta=_META,
                        features=_FS, proximity_bp=200_000).read_text()
    assert "feature notes" in html              # column appears with features
    # chr7 3' end: 45S array 1 kb from the end (700000 - 699000).
    assert "rdna_45S (45S_array_n28) 1 kb from the 3' end" in html


def test_report_states_no_feature_nearby(tmp_path):
    html = write_report(tmp_path / "r.html", _FEAT_CALLS, meta=_META,
                        features=_FS, proximity_bp=200_000).read_text()
    assert "no feature within 200 kb" in html    # chr8, both bare ends


def test_report_capped_end_gets_no_note(tmp_path):
    html = write_report(tmp_path / "r.html", _FEAT_CALLS, meta=_META,
                        features=_FS, proximity_bp=200_000).read_text()
    # chr7's 5' end is capped -> no "5' cap missing" note for it.
    assert "5' cap missing" not in html.split("chr8")[0]
