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
from .prepare import Prepared

# Above this many drawn rectangles we warn (suggest --min-count); still vector.
_DENSE_WARN = 100_000
_BAR_H = 0.3            # bar height in y units (thin backbone; dots sit on top)
_BAR_FILL = "#f0f0f0"  # length-bar body fill in dot style (so the chr is visible)


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
                        labels=None, subtitle=None):
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(labels if labels is not None else order)))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    # Stack title + optional note in *points* above the axes top, so they never
    # collide regardless of figure height (a fig-fraction offset would).
    ax.set_title(f"teloviz — {title}", pad=34 if subtitle else 6)
    if subtitle:
        ax.annotate(subtitle, xy=(0.5, 1), xycoords="axes fraction",
                    xytext=(0, 7), textcoords="offset points",
                    ha="center", va="bottom", fontsize=7.5, color="#555")
    # Margins in *absolute inches* (converted to the fractions subplots_adjust
    # wants), so the title/labels/colorbar keep their real size at ANY figure
    # height — fixed fractions would clip the title on short figures and leave a
    # huge blank on tall ones. The bottom band holds the colorbar.
    fw, fh = fig.get_size_inches()
    left_in, right_in, top_in, bottom_in = 1.0, 0.35, 0.8, 1.35
    fig.subplots_adjust(left=left_in / fw, right=1 - right_in / fw,
                        top=1 - top_in / fh, bottom=bottom_in / fh)
    cax = fig.add_axes((0.34, 0.55 / fh, 0.34, 0.11 / fh))
    cbar = fig.colorbar(scheme.scalar_mappable(), cax=cax, orientation="horizontal",
                        ticks=scheme.tick_positions())
    cbar.ax.set_xticklabels(scheme.tick_labels(), fontsize=7)
    cbar.set_label(scheme.cbar_label())


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
    """
    if call.five:
        ax.scatter([-0.5 * gutter], [y], marker=">", s=34, color="#111", zorder=4)
    if call.three:
        ax.scatter([length + 0.35 * gutter], [y], marker="<", s=34, color="#111", zorder=4)


def render(prepared: Prepared, scheme: ColorScheme, *,
           style: str = "dot", dot_size: float = 40.0, calls=None,
           call_note: str | None = None,
           width: float | None = None, height: float | None = None):
    """Full-length ideogram (one length-proportional bar per chromosome).

    ``width``/``height`` are the figure size in inches (both None → auto: width
    10, height scales with chromosome count). ``calls`` (optional list of
    :class:`~teloviz.calling.Call`) adds telomere-cap markers at capped bar ends
    and a ``*`` on both-ends-capped chromosome labels.
    """
    order, lengths = prepared.order, prepared.lengths
    n = len(order)
    fig_w = width if width else 10.0
    fig_h = height if height else max(2.0, 0.42 * n + 1.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    by_id = {c.id: c for c in calls} if calls else {}
    max_len = max(lengths.values()) if lengths else 1
    # A left gutter (between the labels and x=0) holds the 5' cap triangles, so
    # they never overlap the chromosome names. Only reserved when calls exist.
    gutter = 0.05 * max_len if by_id else 0.0
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
        if cid in by_id:
            _mark_calls(ax, y, lengths[cid], by_id[cid], gutter)
    if style == "dot":
        ax.scatter(dot_xs, dot_ys, c=dot_cs, s=dot_size, edgecolors="none", zorder=3)
    else:
        _draw_rects(ax, rects, facecolors)

    labels = [cid + (" *" if by_id.get(cid) and by_id[cid].both else "") for cid in order]
    ax.set_xlim(-gutter, max_len * 1.01 + 0.5 * gutter)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("Position (Mb)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x / 1e6:g}"))
    _style_and_colorbar(fig, ax, scheme, order, n, scheme.mode,
                        labels=labels, subtitle=call_note)
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
