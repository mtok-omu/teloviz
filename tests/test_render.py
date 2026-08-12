"""Tests for the ideogram renderer (spec sections 6-7).

We drive real ``render`` calls on a tiny prepared table and introspect the
resulting matplotlib Figure: backbone bars, dot offsets / rect paths, axis
labels, and file export. Kept structural (counts / positions), not pixel-exact.
"""

from __future__ import annotations

import matplotlib.colors as mcolors
import pandas as pd
import pytest
from matplotlib.figure import Figure

from teloviz._mpl import plt
from teloviz.color import build_scheme
from teloviz.prepare import prepare
from teloviz.render import render, save

# chr1: window 10000 -> 302 (colored), 20000 -> 0 (white), 30000 -> 3 (colored).
# chr2: window 10000 -> 0 (white), 20000 -> 100 (colored).
# -> 3 non-white sum windows; centres at w - window_size/2 = 5000/25000/15000.
_ROWS = [
    ("chr1", 10000, 300, 2, "T"),
    ("chr1", 20000, 0, 0, "T"),
    ("chr1", 30000, 2, 1, "T"),
    ("chr2", 10000, 0, 0, "T"),
    ("chr2", 20000, 100, 0, "T"),
]
_EXPECTED_DOT_X = [5000.0, 15000.0, 25000.0]  # sorted window centres


def _prepared():
    df = pd.DataFrame(_ROWS, columns=["id", "window", "forward", "reverse", "motif"])
    return prepare(df)


def _scheme(mode="sum"):
    return build_scheme(mode, bin_w=100, cap=500, log=False,
                        cmap_sum="Reds", cmap_div="RdBu_r")


def _main_ax(fig):
    """The ideogram axes (the first one; the colorbar cax is added afterwards)."""
    return fig.axes[0]


@pytest.mark.parametrize("style", ["dot", "rect"])
@pytest.mark.parametrize("mode", ["sum", "orientation"])
def test_render_returns_figure(style, mode):
    fig = render(_prepared(), _scheme(mode), style=style)
    try:
        assert isinstance(fig, Figure)
        # Chromosomes read top -> bottom in natural order; y labels are reversed.
        labels = [t.get_text() for t in _main_ax(fig).get_yticklabels()]
        assert labels == ["chr2", "chr1"]
    finally:
        plt.close(fig)


def test_backbone_one_bar_per_chromosome():
    # Each chromosome contributes exactly one length bar (a Rectangle patch),
    # regardless of style; window marks live in a separate collection.
    for style in ("dot", "rect"):
        fig = render(_prepared(), _scheme("sum"), style=style)
        try:
            assert len(_main_ax(fig).patches) == 2
        finally:
            plt.close(fig)


def test_dot_style_one_point_per_nonwhite_window():
    fig = render(_prepared(), _scheme("sum"), style="dot")
    try:
        offsets = _main_ax(fig).collections[0].get_offsets()
        assert len(offsets) == 3
        xs = sorted(float(x) for x, _ in offsets)
        assert xs == _EXPECTED_DOT_X
    finally:
        plt.close(fig)


def test_dot_size_controls_marker_area():
    fig = render(_prepared(), _scheme("sum"), style="dot", dot_size=123.0)
    try:
        sizes = _main_ax(fig).collections[0].get_sizes()
        assert list(sizes) == [123.0]
    finally:
        plt.close(fig)


def test_rect_style_draws_one_rect_per_nonwhite_window():
    fig = render(_prepared(), _scheme("sum"), style="rect")
    try:
        # ax.collections[0] is the PatchCollection of window rectangles.
        paths = _main_ax(fig).collections[0].get_paths()
        assert len(paths) == 3
    finally:
        plt.close(fig)


def test_calls_asterisk_and_markers():
    from teloviz.calling import Call
    calls = [Call("chr1", 30000, 300, 300, True, True),   # both -> "chr1 *"
             Call("chr2", 20000, 300, 0, True, False)]    # 5' only -> no asterisk
    prepared, scheme = _prepared(), _scheme("sum")
    plain = render(prepared, scheme, style="dot")
    marked = render(prepared, scheme, style="dot", calls=calls,
                    call_note="telomere call")
    try:
        labels = [t.get_text() for t in _main_ax(marked).get_yticklabels()]
        assert "chr1 *" in labels and "chr2" in labels and "chr2 *" not in labels
        # The capped-end triangles are extra collections beyond the dot scatter.
        assert len(_main_ax(marked).collections) > len(_main_ax(plain).collections)
    finally:
        plt.close(plain); plt.close(marked)


def _amber_marks(ax):
    """Collections drawn in the warning amber (one per flagged capped end)."""
    warn = mcolors.to_rgba("#b8860b")
    return [c for c in ax.collections
            if len(c.get_facecolor()) and tuple(c.get_facecolor()[0]) == warn]


def test_short_contig_call_is_ambered_and_parenthesized():
    from teloviz.calling import Call
    # chr1 flagged short (end regions overlapped), chr2 a normal both-ends call.
    calls = [Call("chr1", 30000, 300, 300, True, True, short=True),
             Call("chr2", 20000, 300, 300, True, True)]
    fig = render(_prepared(), _scheme("sum"), style="dot", calls=calls)
    try:
        ax = _main_ax(fig)
        ticks = {t.get_text(): t for t in ax.get_yticklabels()}
        assert set(ticks) == {"chr1 (*)", "chr2 *"}   # uncertain call parenthesized
        warn = mcolors.to_rgba("#b8860b")
        assert mcolors.to_rgba(ticks["chr1 (*)"].get_color()) == warn
        assert mcolors.to_rgba(ticks["chr2 *"].get_color()) != warn
        # Exactly two amber triangles (chr1's two ends); chr2's stay black.
        assert len(_amber_marks(ax)) == 2
    finally:
        plt.close(fig)


def test_uncalled_short_contig_is_not_flagged():
    # Short but nothing called: no verdict to distrust, so the name stays black.
    from teloviz.calling import Call
    calls = [Call("chr1", 30000, 0, 0, False, False, short=True),
             Call("chr2", 20000, 300, 0, True, False)]
    fig = render(_prepared(), _scheme("sum"), style="dot", calls=calls)
    try:
        ticks = {t.get_text(): t for t in _main_ax(fig).get_yticklabels()}
        assert mcolors.to_rgba(ticks["chr1"].get_color()) != mcolors.to_rgba("#b8860b")
    finally:
        plt.close(fig)


def _feature_set():
    from teloviz.features import Feature, FeatureSet
    # chr1 gets one rdna_45S feature; chr2 has none.
    return FeatureSet(
        features=[Feature("chr1", 5000, 8000, "arr", "rdna_45S")],
        lane_order=["rdna_45S"],
    )


def test_features_add_lane_patches_below_the_bar():
    from teloviz.features import color_for_type
    prepared, scheme = _prepared(), _scheme("sum")
    plain = render(prepared, scheme, style="dot")
    withf = render(prepared, scheme, style="dot", features=_feature_set())
    try:
        # Exactly one extra patch (the feature rect) beyond the 2 backbone bars.
        assert len(_main_ax(withf).patches) == len(_main_ax(plain).patches) + 1
        green = color_for_type("rdna_45S", ["rdna_45S"])
        feats = [p for p in _main_ax(withf).patches
                 if mcolors.to_hex(p.get_facecolor()[:3]) == green.lower()]
        assert len(feats) == 1
        # The feature sits strictly below the bar it annotates. chr1 is the top
        # row (order [chr1, chr2] -> y=1); the bar spans y +/- 0.15, so the lane
        # top must be at or below y - 0.15 = 0.85.
        assert feats[0].get_y() + feats[0].get_height() <= 1 - 0.15 + 1e-9
    finally:
        plt.close(plain); plt.close(withf)


def test_features_do_not_change_bar_colors():
    # The window dot colors must be identical with and without features.
    prepared, scheme = _prepared(), _scheme("sum")
    a = render(prepared, scheme, style="dot")
    b = render(prepared, scheme, style="dot", features=_feature_set())
    try:
        ca = _main_ax(a).collections[0].get_facecolor()
        cb = _main_ax(b).collections[0].get_facecolor()
        assert (ca == cb).all()
    finally:
        plt.close(a); plt.close(b)


def test_tiny_feature_widened_and_right_clamped():
    # A 1 bp feature flush against chr2's right end (length 20000) must be widened
    # to the min visible width yet never overflow past the chromosome end.
    from teloviz.features import Feature, FeatureSet
    fs = FeatureSet(features=[Feature("chr2", 19999, 20000, "tip", "rdna_45S")],
                    lane_order=["rdna_45S"])
    fig = render(_prepared(), _scheme("sum"), style="dot", features=fs)
    try:
        feats = [p for p in _main_ax(fig).patches
                 if abs(p.get_height()) < 0.2 and p.get_width() > 1.5]
        assert feats, "tiny feature should be widened to be visible"
        f = feats[0]
        assert f.get_x() + f.get_width() <= 20000 + 1e-6   # clamped, no overflow
    finally:
        plt.close(fig)


def test_save_writes_every_requested_format(tmp_path):
    fig = render(_prepared(), _scheme("sum"), style="dot")
    try:
        prefix = str(tmp_path / "out")
        paths = save(fig, prefix, "sum", ["pdf", "png"], dpi=100)
    finally:
        plt.close(fig)
    assert [p.name for p in paths] == ["out.sum.pdf", "out.sum.png"]
    for p in paths:
        assert p.exists() and p.stat().st_size > 0
