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


def _style_and_colorbar(fig, ax, scheme: ColorScheme, order, n, title):
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(order)))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title(f"teloviz — {title}")
    # A dedicated bottom colorbar axes (figure coords) keeps a fixed, small gap
    # regardless of how tall the ideogram is — fig.colorbar(ax=...) would reserve
    # a fraction of a very tall axes and leave a huge blank.
    fig.subplots_adjust(left=0.12, right=0.97, top=0.95, bottom=0.09)
    cax = fig.add_axes((0.30, 0.045, 0.40, 0.010))
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


def render(prepared: Prepared, scheme: ColorScheme, *,
           style: str = "dot", dot_size: float = 40.0,
           width: int | None = None, height: int | None = None):
    """Full-length ideogram (one length-proportional bar per chromosome)."""
    order, lengths = prepared.order, prepared.lengths
    n = len(order)
    fig_w = width / 100 if width else 10.0
    fig_h = height / 100 if height else max(2.0, 0.42 * n + 1.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    max_len = max(lengths.values()) if lengths else 1
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
    if style == "dot":
        ax.scatter(dot_xs, dot_ys, c=dot_cs, s=dot_size, edgecolors="none", zorder=3)
    else:
        _draw_rects(ax, rects, facecolors)

    ax.set_xlim(0, max_len * 1.01)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("Position (Mb)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x / 1e6:g}"))
    _style_and_colorbar(fig, ax, scheme, order, n, scheme.mode)
    return fig


def save(fig, out_prefix: str, label: str, formats: list[str], dpi: int) -> list[Path]:
    """Export the figure to every requested format; return the paths written."""
    paths: list[Path] = []
    for fmt in formats:
        p = Path(f"{out_prefix}.{label}.{fmt}")
        if p.parent != Path(""):
            p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        paths.append(p)
    return paths
