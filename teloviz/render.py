"""Ideogram rendering + export for teloviz (spec sections 6-7).

Each chromosome is one length-proportional horizontal strip. A strip is drawn as
a single ``imshow`` row (one column per window) rather than thousands of patches,
so a 20-chromosome / 200k-window genome renders fast and stays small. Missing
windows and the unmeasured tail are NaN -> white. A thin outline marks the true
sequence length; data-free regions stay white.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from matplotlib.patches import Rectangle
from matplotlib.ticker import FuncFormatter

from ._mpl import plt
from .color import ColorScheme
from .prepare import Prepared


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

    for i, cid in enumerate(order):
        y = n - 1 - i  # first (natural-sorted) id on top
        length = lengths[cid]
        ncols = max(1, int(np.ceil(length / ws)))

        grid = np.full(ncols, np.nan)
        sub = table[table["id"] == cid]
        col = np.rint(sub["window"].to_numpy() / ws).astype(int) - 1
        vals = scheme.values(sub["forward"].to_numpy(), sub["reverse"].to_numpy())
        valid = (col >= 0) & (col < ncols)
        grid[col[valid]] = scheme.transform(vals[valid])

        data = np.ma.masked_invalid(grid).reshape(1, -1)
        ax.imshow(
            data, aspect="auto", cmap=scheme.cmap, norm=scheme.norm,
            extent=[0, ncols * ws, y - 0.4, y + 0.4],
            interpolation="nearest", origin="lower", zorder=1,
        )
        # True-length outline so white regions never hide the bar's extent/ends.
        ax.add_patch(Rectangle((0, y - 0.4), length, 0.8, fill=False,
                               edgecolor="black", linewidth=0.5, zorder=2))

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
