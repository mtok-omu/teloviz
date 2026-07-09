"""Ideogram rendering + export for teloviz (spec sections 6-7).

Each chromosome is one length-proportional horizontal bar. Only *non-white*
windows (telomere arrays, internal clusters) are drawn, as true-position vector
rectangles — so the bar stays honestly length-proportional (no fattening) and
the PDF/SVG is zoomable in any viewer: zoom into an end and the thin telomere
rectangle stays crisp at its real size. Background/missing windows are simply not
drawn (white). A thin outline marks the true sequence length.
"""

from __future__ import annotations

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


def render(
    prepared: Prepared,
    scheme: ColorScheme,
    *,
    width: int | None = None,
    height: int | None = None,
):
    """Build (but do not save) a matplotlib Figure for one mode."""
    order = prepared.order
    lengths = prepared.lengths
    ws = prepared.window_size
    n = len(order)

    fig_w = width / 100 if width else 10.0
    fig_h = height / 100 if height else max(2.0, 0.42 * n + 1.6)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))

    max_len = max(lengths.values()) if lengths else 1
    table = prepared.table

    rects: list[Rectangle] = []
    facecolors: list = []
    for i, cid in enumerate(order):
        y = n - 1 - i  # first (natural-sorted) id on top
        length = lengths[cid]

        sub = table[table["id"] == cid]
        windows = sub["window"].to_numpy()
        vals = scheme.values(sub["forward"].to_numpy(), sub["reverse"].to_numpy())
        colors = scheme.rgba(vals)
        for w, rgba in zip(windows, colors):
            if np.allclose(rgba, (1.0, 1.0, 1.0, 1.0)):
                continue  # background / net~0 / below --min-count -> draw nothing
            start = max(0.0, w - ws)
            end = min(float(w), length)
            if end <= start:
                continue
            rects.append(Rectangle((start, y - 0.4), end - start, 0.8))
            facecolors.append(rgba)

        # True-length outline so the bar's extent and ends read even when white.
        ax.add_patch(Rectangle((0, y - 0.4), length, 0.8, fill=False,
                               edgecolor="black", linewidth=0.5, zorder=2))

    if len(rects) > _DENSE_WARN:
        import sys
        print(
            f"teloviz: warning: drawing {len(rects)} window rectangles; consider "
            f"--min-count to suppress background and shrink the figure.",
            file=sys.stderr,
        )
    if rects:
        ax.add_collection(PatchCollection(rects, facecolors=facecolors,
                                          edgecolors="none", zorder=1))

    ax.set_xlim(0, max_len * 1.01)
    ax.set_ylim(-0.6, n - 0.4)
    ax.set_yticks(range(n))
    ax.set_yticklabels(list(reversed(order)))
    ax.set_xlabel("Position (Mb)")
    ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _pos: f"{x / 1e6:g}"))
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.set_title(f"teloviz — {scheme.mode}")

    cbar = fig.colorbar(
        scheme.scalar_mappable(), ax=ax, orientation="horizontal",
        fraction=0.06, pad=0.32, aspect=40, ticks=scheme.tick_positions(),
    )
    cbar.ax.set_xticklabels(scheme.tick_labels(), fontsize=7)
    cbar.set_label(scheme.cbar_label())

    fig.tight_layout()
    return fig


def save(fig, out_prefix: str, mode: str, formats: list[str], dpi: int) -> list[Path]:
    """Export the figure to every requested format; return the paths written."""
    paths: list[Path] = []
    for fmt in formats:
        p = Path(f"{out_prefix}.{mode}.{fmt}")
        if p.parent != Path(""):
            p.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(p, dpi=dpi, bbox_inches="tight")
        paths.append(p)
    return paths
