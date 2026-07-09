"""Ideogram rendering + export for teloviz (spec sections 6-7).

Two layouts, both drawing only *non-white* windows (telomere arrays / internal
clusters) as true-position vector rectangles — honest length-proportional, no
fattening, and zoomable in any PDF/SVG viewer:

- ``render``      : full-length bars (best for internal clusters / overview).
- ``render_ends`` : a compact both-ends view. Each chromosome shows its first N
  and last N bp joined by a ``≈`` break, so the elided middle costs no width.
  Good for publication and for reading which ends are telomere-capped (T2T).
"""

from __future__ import annotations

import math
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
_BAR_H = 0.8  # bar height in y units


def _nice_step(span: float) -> float:
    """A round tick step giving ~3-4 ticks across ``span``."""
    raw = span / 3 if span > 0 else 1
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 2.5, 5, 10):
        if m * mag >= raw:
            return m * mag
    return 10 * mag


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


def render(prepared: Prepared, scheme: ColorScheme, *,
           width: int | None = None, height: int | None = None):
    """Full-length ideogram (one length-proportional bar per chromosome)."""
    order, lengths = prepared.order, prepared.lengths
    n = len(order)
    fig_w = width / 100 if width else 10.0
    fig_h = height / 100 if height else max(2.0, 0.42 * n + 1.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    max_len = max(lengths.values()) if lengths else 1
    rects, facecolors = [], []
    for i, cid in enumerate(order):
        y = n - 1 - i
        for start, end, rgba in _iter_windows(prepared, scheme, cid):
            if end > start:
                rects.append(Rectangle((start, y - _BAR_H / 2), end - start, _BAR_H))
                facecolors.append(rgba)
        ax.add_patch(Rectangle((0, y - _BAR_H / 2), lengths[cid], _BAR_H,
                               fill=False, edgecolor="black", linewidth=0.5, zorder=2))
    _draw_rects(ax, rects, facecolors)

    ax.set_xlim(0, max_len * 1.01)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_xlabel("Position (Mb)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _p: f"{x / 1e6:g}"))
    _style_and_colorbar(fig, ax, scheme, order, n, scheme.mode)
    return fig


def render_ends(prepared: Prepared, scheme: ColorScheme, *, ends_bp: int,
                width: int | None = None, height: int | None = None):
    """Both-ends view: first/last ``ends_bp`` per chromosome, joined by a break."""
    order, lengths = prepared.order, prepared.lengths
    n = len(order)
    N = float(ends_bp)
    gap = 0.10 * N          # visual break width
    offset = N + gap        # x where the right (3') panel starts
    total_x = 2 * N + gap

    fig_w = width / 100 if width else 8.0
    fig_h = height / 100 if height else max(2.0, 0.42 * n + 1.8)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    rects, facecolors = [], []
    for i, cid in enumerate(order):
        y = n - 1 - i
        L = lengths[cid]
        y0 = y - _BAR_H / 2
        for start, end, rgba in _iter_windows(prepared, scheme, cid):
            # Left (5') panel: genomic [0, N] -> x = genomic.
            ls, le = start, min(end, N)
            if le > ls and ls < N:
                rects.append(Rectangle((ls, y0), le - ls, _BAR_H)); facecolors.append(rgba)
            # Right (3') panel: genomic [L-N, L] -> x = offset + (g - (L-N)).
            rs, re = max(start, L - N), end
            if re > rs and re > L - N:
                x0 = offset + (rs - (L - N))
                rects.append(Rectangle((x0, y0), re - rs, _BAR_H)); facecolors.append(rgba)
        # Panel outlines (only the part that exists for short chromosomes).
        ax.add_patch(Rectangle((0, y0), min(N, L), _BAR_H, fill=False,
                               edgecolor="black", linewidth=0.5, zorder=2))
        ax.add_patch(Rectangle((offset + max(0.0, N - L), y0), min(N, L), _BAR_H,
                               fill=False, edgecolor="black", linewidth=0.5, zorder=2))
        # Break marker in the gap.
        if L > 2 * N:
            ax.text(N + gap / 2, y, "≈", ha="center", va="center", fontsize=9)
    _draw_rects(ax, rects, facecolors)

    ax.set_xlim(-0.02 * N, total_x + 0.02 * N)
    ax.set_ylim(-0.6, n - 0.4)

    step = _nice_step(N)
    left_pos = np.arange(0, N + 1, step)
    right_pos = [offset + (N - d) for d in left_pos]          # d = Mb from 3' end
    ax.set_xticks(list(left_pos) + list(right_pos))
    ax.set_xticklabels([f"{p / 1e6:g}" for p in left_pos]
                       + [f"{d / 1e6:g}" for d in left_pos], fontsize=8)
    # Panel captions, placed a fixed offset below the axis (not axes-fraction, so
    # they stay put on tall figures).
    for x, cap in ((N / 2, "Mb from 5′ end"), (offset + N / 2, "Mb from 3′ end")):
        ax.annotate(cap, xy=(x, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, -24), textcoords="offset points",
                    ha="center", va="top")

    _style_and_colorbar(fig, ax, scheme, order, n, f"{scheme.mode} · ends {N / 1e6:g} Mb")
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
