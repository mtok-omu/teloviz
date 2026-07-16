"""Tests for the optional annotation-track features (annotation-track spec)."""

from __future__ import annotations

from teloviz.features import (
    color_for_type,
    lane_label,
    load_features,
    normalize_features,
)

_BED = (
    "# a comment line\n"
    "\n"
    "chr7\t690000\t699000\t45S_array_n28\trdna_45S\n"
    "chr7\t300000\t318000\t5S_array_n42\trdna_5S\n"
    "chr1\t100\t200\n"                       # cols 4/5 default
)


def _write(tmp_path, text):
    p = tmp_path / "features.bed"
    p.write_text(text)
    return p


def test_load_parses_columns_and_defaults(tmp_path):
    feats = load_features(_write(tmp_path, _BED))
    assert len(feats) == 3
    f0 = feats[0]
    assert (f0.chrom, f0.start, f0.end, f0.name, f0.type) == (
        "chr7", 690000, 699000, "45S_array_n28", "rdna_45S")
    # Missing name/type columns fall back to defaults.
    f2 = feats[2]
    assert f2.name == "chr1:100-200" and f2.type == "feature"


def test_load_rejects_too_few_columns(tmp_path):
    import pytest
    with pytest.raises(ValueError, match=">= 3"):
        load_features(_write(tmp_path, "chr1\t100\n"))


def test_normalize_drops_unknown_chrom(tmp_path, capsys):
    feats = load_features(_write(tmp_path, _BED + "chrZ\t1\t2\tx\trdna_5S\n"))
    fs = normalize_features(feats, {"chr7": 700000, "chr1": 1000})
    assert all(f.chrom in ("chr7", "chr1") for f in fs.features)
    assert "dropped 1 feature" in capsys.readouterr().err


def test_normalize_drops_empty_interval(tmp_path):
    feats = load_features(_write(tmp_path, "chr1\t500\t400\tbad\trdna_5S\n"))
    fs = normalize_features(feats, {"chr1": 1000}, quiet=True)
    assert fs.features == []


def test_normalize_clamps_end_to_length(tmp_path):
    feats = load_features(_write(tmp_path, "chr1\t900\t5000\tspill\trdna_5S\n"))
    fs = normalize_features(feats, {"chr1": 1000}, quiet=True)
    assert len(fs.features) == 1 and fs.features[0].end == 1000


def test_lane_order_is_first_appearance(tmp_path):
    feats = load_features(_write(tmp_path, _BED))
    fs = normalize_features(feats, {"chr7": 700000, "chr1": 1000}, quiet=True)
    # rdna_45S appears before rdna_5S, then the defaulted "feature".
    assert fs.lane_order == ["rdna_45S", "rdna_5S", "feature"]


def test_color_map_and_palette_fallback():
    order = ["rdna_45S", "mynovel"]
    assert color_for_type("rdna_45S", order) == "#2E7D32"
    # Unknown type gets a palette color (stable, not a mapped one).
    assert color_for_type("mynovel", order).startswith("#")
    assert color_for_type("mynovel", order) != "#2E7D32"


def test_candidate_label_suffix():
    assert lane_label("rdna_45S") == "rdna_45S"
    assert lane_label("centromere") == "centromere (candidate)"
