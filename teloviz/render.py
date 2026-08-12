"""Ideogram rendering + export for teloviz (spec sections 6-7).

``render`` draws a full-length ideogram: one length-proportional bar per
chromosome, showing only *non-white* windows (telomere arrays / internal
clusters = misjoin hints).

Two mark styles (``style=``):

- ``dot``  (default): each non-white window is a fixed-size round marker at its
  genomic centre on a light-gray length bar. Telomere arrays sit in a handful of
  tiny windows, so honest-width rectangles vanish to invisible slivers on a
  multi-Mb bar; a fixed on-screen dot size (``--dot-size``) stays readable
  regardless of window width. Position stays true; only the mark size is fixed.
- ``rect``: the honest length-proportional filled rectangle (true width, no
  fattening, zoomable in any PDF/SVG viewer). Kept for internal-cluster extent.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

from ._mpl import plt
from .color import ColorScheme
from .features import FEATURE_ALPHA, FeatureSet, color_for_type, lane_label
from .prepare import Prepared

# Above this many drawn rectangles we warn (suggest --min-count); still vector.
_DENSE_WARN = 100_000
_BAR_H = 0.3            # bar height in y units (thin backbone; dots sit on top)
_BAR_FILL = "#f0f0f0"  # length-bar body fill in dot style (so the chr is visible)
_CALL_COLOR = "#111"   # telomere-cap markers on a chromosome-length sequence
# Sequences too short for the two end regions to be disjoint get their call marks
# in this amber (the same warning color the HTML report uses) instead of black.
_WARN_COLOR = "#b8860b"

# Feature lane geometry, in the same y units as the bars (bars are 1.0 apart).
# All feature types share ONE thin lane just under the bar, distinguished by
# color, so a feature never drifts down toward the next chromosome (which made it
# ambiguous which bar it belonged to). The lane hugs the bar it annotates.
_LANE_TOP_GAP = 0.05   # bar bottom -> lane top
_LANE_H = 0.14         # lane height


def _iter_windows(prepared: Prepared, scheme: ColorScheme, cid: str):
    """Yield (genomic_start, genomic_end, rgba) for non-white windows of ``cid``."""
    ws = prepared.window_size
    length = prepared.lengths[cid]
    sub = prepared.table[prepared.table["id"] == cid]
    vals = scheme.values(sub["forward"].to_numpy(), sub["reverse"].to_numpy())
    colors = scheme.rgba(vals)
    for w, rgba in zip(sub["window"].to_numpy(), colors):
        if np.allclose(rgba, (1.0, 1.0, 1.0, 1.0)):
            continue
        yield max(0.0, w - ws), min(float(w), length), rgba


def _draw_rects(ax, rects, facecolors):
    if len(rects) > _DENSE_WARN:
        print(
            f"teloviz: warning: drawing {len(rects)} window rectangles; consider "
            f"--min-count to suppress background and shrink the figure.",
            file=sys.stderr,
        )
    if rects:
        ax.add_collection(PatchCollection(rects, facecolors=facecolors,
                                          edgecolors="none", zorder=1))


def _style_and_colorbar(fig, ax, scheme: ColorScheme, order, n, title,
                        labels=None, subtitle=None, fs=10.0, feat_legend=None):
    # All text sizes are set explicitly (as multiples of the base font size ``fs``)
    # so they survive save() regardless of rcParams. fs=10 reproduces the defaults.
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(labels if labels is not None else order)),
                       fontsize=fs)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    # Stack title + optional note in *points* above the axes top, so they never
    # collide regardless of figure height (a fig-fraction offset would).
    ax.set_title(f"teloviz — {title}", pad=3.4 * fs if subtitle else 6,
                 fontsize=1.2 * fs)
    if subtitle:
        ax.annotate(subtitle, xy=(0.5, 1), xycoords="axes fraction",
                    xytext=(0, 0.7 * fs), textcoords="offset points",
                    ha="center", va="bottom", fontsize=0.75 * fs, color="#555")
    # Margins in *absolute inches* (converted to the fractions subplots_adjust
    # wants), so the title/labels/colorbar keep their real size at ANY figure
    # height — fixed fractions would clip the title on short figures and leave a
    # huge blank on tall ones. Margins scale with the font so larger text still
    # fits. The bottom band holds the colorbar.
    fw, fh = fig.get_size_inches()
    r = fs / 10.0
    left_in, right_in, top_in, bottom_in = 1.0 * r, 0.35, 0.8 * r, 1.35 * r
    # Clamp so the margins can never exceed the figure (a small figure with a big
    # font would otherwise make bottom >= top and crash); shrink them to fit.
    if top_in + bottom_in > 0.85 * fh:
        s = 0.85 * fh / (top_in + bottom_in)
        top_in *= s; bottom_in *= s
    if left_in + right_in > 0.85 * fw:
        s = 0.85 * fw / (left_in + right_in)
        left_in *= s; right_in *= s
    fig.subplots_adjust(left=left_in / fw, right=1 - right_in / fw,
                        top=1 - top_in / fh, bottom=bottom_in / fh)
    cax = fig.add_axes((0.34, 0.55 * r / fh, 0.34, 0.11 / fh))
    cbar = fig.colorbar(scheme.scalar_mappable(), cax=cax, orientation="horizontal",
                        ticks=scheme.tick_positions())
    cbar.ax.set_xticklabels(scheme.tick_labels(), fontsize=0.7 * fs)
    cbar.set_label(scheme.cbar_label(), fontsize=fs)
    # Feature-lane key in the bottom-left margin (left of the colorbar), one swatch
    # per type. A single legend replaces per-row lane labels, which would be
    # unreadable repeated on every chromosome in a static multi-chrom ideogram.
    if feat_legend:
        from matplotlib.patches import Patch
        handles = [Patch(facecolor=c, edgecolor="none", alpha=0.6, label=lbl)
                   for lbl, c in feat_legend]
        fig.legend(handles=handles, loc="lower left",
                   bbox_to_anchor=(left_in / fw, 0.3 * r / fh),
                   frameon=False, fontsize=0.75 * fs, title="features",
                   title_fontsize=0.8 * fs, handlelength=1.2, handleheight=1.0,
                   borderaxespad=0.0)


def _add_backbone(ax, x0, y, width, style):
    """Draw one chromosome length bar; filled gray for dots, outline for rects."""
    if style == "dot":
        ax.add_patch(Rectangle((x0, y - _BAR_H / 2), width, _BAR_H,
                               facecolor=_BAR_FILL, edgecolor="black",
                               linewidth=0.5, zorder=1))
    else:
        ax.add_patch(Rectangle((x0, y - _BAR_H / 2), width, _BAR_H, fill=False,
                               edgecolor="black", linewidth=0.5, zorder=2))


def _mark_calls(ax, y, length, call, gutter):
    """Inward triangles marking telomere-capped ends of one chromosome bar.

    The 5' triangle sits in the left *gutter* — the strip between the y-axis
    (chromosome labels) and x=0 where the bars start — so it can never touch the
    label. The 3' triangle sits just past the bar's right end, clear of the bar
    and its dots. Both stay inside the axes (no clipping into the label margin).

    On a sequence too short for the end regions to be disjoint (``call.short``)
    the triangles are drawn amber, because the same windows fed both ends and the
    call may be spurious.
    """
    color = _WARN_COLOR if call.short else _CALL_COLOR
    if call.five:
        ax.scatter([-0.5 * gutter], [y], marker=">", s=34, color=color, zorder=4)
    if call.three:
        ax.scatter([length + 0.35 * gutter], [y], marker="<", s=34, color=color, zorder=4)


def _label_suffix(call) -> str:
    """``Call.suffix`` for a called chromosome, nothing for an uncalled one."""
    return call.suffix if call is not None else ""


def _draw_feature_lanes(ax, y, cid, features: FeatureSet, length, min_w):
    """Draw this chromosome's features as neutral rects in one shared lane just
    below the bar, colored by ``type``. Bar colors are never touched; features
    live strictly outside the bar. Tiny features are widened to ``min_w`` and
    clamped to the chromosome ends so a 9.8 kb array is still visible on a 67 Mb
    bar and never overflows the right end. Wider features are drawn first so a
    narrow feature of another type stacked at the same locus stays on top.
    """
    top = y - _BAR_H / 2 - _LANE_TOP_GAP
    feats = sorted(features.for_chrom(cid), key=lambda f: f.end - f.start,
                   reverse=True)
    for f in feats:
        w = max(f.end - f.start, min_w)
        x0 = f.start
        if x0 + w > length:      # right-end clamp: extend leftward, stay flush
            x0 = length - w
        if x0 < 0:               # feature wider than the whole (tiny) chromosome
            x0, w = 0.0, min(w, length)
        color = color_for_type(f.type, features.lane_order)
        ax.add_patch(Rectangle((x0, top - _LANE_H), w, _LANE_H, facecolor=color,
                               edgecolor="none", alpha=FEATURE_ALPHA, zorder=2))


def render(prepared: Prepared, scheme: ColorScheme, *,
           style: str = "dot", dot_size: float = 40.0, calls=None,
           call_note: str | None = None, font_size: float = 10.0,
           width: float | None = None, height: float | None = None,
           features: FeatureSet | None = None):
    """Full-length ideogram (one length-proportional bar per chromosome).

    ``width``/``height`` are the figure size in inches (both None → auto: width
    10, height scales with chromosome count). ``font_size`` is the base text size
    in points (all text scales from it; 10 reproduces the defaults). ``calls``
    (optional list of :class:`~teloviz.calling.Call`) adds telomere-cap markers at
    capped bar ends and a ``*`` on both-ends-capped chromosome labels. A call on a
    sequence flagged ``short`` (end regions overlap) is drawn amber and its ``*``
    parenthesized, since it may be spurious.
    """
    order, lengths = prepared.order, prepared.lengths
    n = len(order)
    r = font_size / 10.0
    fig_w = width if width else 10.0
    # Auto height scales with the font so big labels get taller rows/margins.
    fig_h = height if height else max(2.5 * r, (0.42 * n + 1.6) * r)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    by_id = {c.id: c for c in calls} if calls else {}
    max_len = max(lengths.values()) if lengths else 1
    # A left gutter (between the labels and x=0) holds the 5' cap triangles, so
    # they never overlap the chromosome names. Only reserved when calls exist.
    gutter = 0.05 * max_len if by_id else 0.0
    # Feature track: one shared lane just below each bar, colored by type. Tiny
    # features are widened to a fraction of the longest bar so they stay visible;
    # positions/right-edge stay honest via clamping.
    has_feat = bool(features and features.lane_order)
    min_feat_w = max_len / 400.0
    rects, facecolors = [], []
    dot_xs, dot_ys, dot_cs = [], [], []
    for i, cid in enumerate(order):
        y = n - 1 - i
        _add_backbone(ax, 0, y, lengths[cid], style)
        for start, end, rgba in _iter_windows(prepared, scheme, cid):
            if style == "dot":
                dot_xs.append((start + end) / 2.0); dot_ys.append(y); dot_cs.append(rgba)
            elif end > start:
                rects.append(Rectangle((start, y - _BAR_H / 2), end - start, _BAR_H))
                facecolors.append(rgba)
        if has_feat:
            _draw_feature_lanes(ax, y, cid, features, lengths[cid], min_feat_w)
        if cid in by_id:
            _mark_calls(ax, y, lengths[cid], by_id[cid], gutter)
    if style == "dot":
        ax.scatter(dot_xs, dot_ys, c=dot_cs, s=dot_size, edgecolors="none", zorder=3)
    else:
        _draw_rects(ax, rects, facecolors)

    labels = [cid + _label_suffix(by_id.get(cid)) for cid in order]
    ax.set_xlim(-gutter, max_len * 1.01 + 0.5 * gutter)
    # Extend the bottom so the lowest chromosome's feature lane is never clipped.
    stack = (_LANE_TOP_GAP + _LANE_H) if has_feat else 0.0
    ax.set_ylim(min(-0.6, -(_BAR_H / 2 + stack + 0.08)), n - 0.4)
    ax.set_xlabel("Position (Mb)", fontsize=font_size)
    ax.tick_params(axis="x", labelsize=font_size)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x / 1e6:g}"))
    feat_legend = ([(lane_label(t), color_for_type(t, features.lane_order))
                    for t in features.lane_order] if has_feat else None)
    _style_and_colorbar(fig, ax, scheme, order, n, scheme.mode,
                        labels=labels, subtitle=call_note, fs=font_size,
                        feat_legend=feat_legend)
    # Amber the names of flagged sequences so the warning is visible even where
    # the marks are (a single triangle in the gutter is easy to miss). Only rows
    # that actually carry a call are flagged: with nothing called there is no
    # verdict to distrust. Tick labels run bottom -> top, i.e. reversed(order).
    warn_ids = {c.id for c in by_id.values() if c.short and (c.five or c.three)}
    if warn_ids:
        for tick, cid in zip(ax.get_yticklabels(), reversed(order)):
            if cid in warn_ids:
                tick.set_color(_WARN_COLOR)
    return fig


def save(fig, out_prefix: str, label: str, formats: list[str], dpi: int,
         tight: bool = True) -> list[Path]:
    """Export the figure to every requested format; return the paths written.

    ``tight`` trims surrounding whitespace (nice by default) but then the saved
    aspect ratio follows the content, not the requested figure size. Pass
    ``tight=False`` to honor the exact figure dimensions (used when the user sets
    both --width and --height, so the requested ratio is preserved precisely).
    """
    bbox = "tight" if tight else None
    paths: list[Path] = []
    for fmt in formats:
        p = Path(f"{out_prefix}.{label}.{fmt}")
        if p.parent != Path(""):
            p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=dpi, bbox_inches=bbox)
        paths.append(p)
    return paths
