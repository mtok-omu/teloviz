"""Value -> color mapping for teloviz (spec section 5).

Two modes, both discretized into ``--bin``-wide bands and clamped at ``--cap``:

- sum:         v = forward + reverse, white -> red (sequential). v == 0 is its
               own white band, so only telomeric regions light up.
- orientation: v = forward - reverse (net), blue <- white -> red (diverging),
               centered at 0.

``--log`` maps values in log space (log1p for sum; sign * log1p(|v|) for
orientation). We build the *same* discrete cmap + BoundaryNorm used both for the
window colors and the labeled colorbar, so they always agree.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import BoundaryNorm, ListedColormap

from ._mpl import plt


def _sum_edges(bin_w: int, cap: float) -> np.ndarray:
    """Band edges: [0, 1, bin, 2*bin, ..., cap]. The [0,1) band is exact zero."""
    edges = [0.0, 1.0]
    x = float(bin_w)
    while x < cap:
        edges.append(x)
        x += bin_w
    edges.append(float(cap))
    return np.array(sorted(set(edges)), dtype=float)


def _orient_edges(bin_w: int, cap: float) -> np.ndarray:
    """Symmetric band edges: [-cap, ..., -bin, 0, bin, ..., cap]."""
    pos = list(np.arange(bin_w, cap, bin_w)) + [float(cap)]
    pos = sorted(set(float(p) for p in pos if p <= cap))
    return np.array([-p for p in reversed(pos)] + [0.0] + pos, dtype=float)


@dataclass
class ColorScheme:
    mode: str
    edges: np.ndarray            # original-space band edges
    edges_t: np.ndarray          # transformed-space band edges (== edges if linear)
    cmap: ListedColormap
    norm: BoundaryNorm
    log: bool
    cap: float

    def values(self, forward: np.ndarray, reverse: np.ndarray) -> np.ndarray:
        forward = np.asarray(forward, dtype=float)
        reverse = np.asarray(reverse, dtype=float)
        return forward + reverse if self.mode == "sum" else forward - reverse

    def transform(self, v: np.ndarray) -> np.ndarray:
        """Map raw values into the transformed space the norm expects."""
        v = np.asarray(v, dtype=float)
        if self.mode == "sum":
            v = np.clip(v, 0.0, self.cap)
            return np.log1p(v) if self.log else v
        v = np.clip(v, -self.cap, self.cap)
        return np.sign(v) * np.log1p(np.abs(v)) if self.log else v

    def scalar_mappable(self) -> ScalarMappable:
        sm = ScalarMappable(norm=self.norm, cmap=self.cmap)
        sm.set_array([])
        return sm

    def tick_positions(self) -> np.ndarray:
        return self.edges_t

    def tick_labels(self) -> list[str]:
        return [f"{e:g}" for e in self.edges]

    def cbar_label(self) -> str:
        base = "forward + reverse" if self.mode == "sum" else "forward - reverse (net)"
        return f"{base}{' (log)' if self.log else ''}"


def build_scheme(
    mode: str,
    *,
    bin_w: int,
    cap: float,
    log: bool,
    cmap_sum: str,
    cmap_div: str,
) -> ColorScheme:
    if mode == "sum":
        edges = _sum_edges(bin_w, cap)
        n_bands = len(edges) - 1
        base = plt.get_cmap(cmap_sum)
        # Band 0 (exact zero) is white; the rest ramp up the sequential map.
        colors = np.vstack([[1, 1, 1, 1], base(np.linspace(0.15, 1.0, n_bands - 1))])
    else:
        edges = _orient_edges(bin_w, cap)
        n_bands = len(edges) - 1
        base = plt.get_cmap(cmap_div)
        colors = base(np.linspace(0.0, 1.0, n_bands))
        # Force the two bands straddling 0 (|net| < bin) to pure white, so the
        # chromosome body reads as neutral and only clear strand bias colors.
        zero = int(np.argmin(np.abs(edges)))  # index of the 0.0 edge
        for b in (zero - 1, zero):
            if 0 <= b < n_bands:
                colors[b] = [1, 1, 1, 1]

    # NaN (unmeasured tail / gaps) -> white.
    cmap = ListedColormap(colors).with_extremes(bad="white")

    if log:
        if mode == "sum":
            edges_t = np.log1p(edges)
        else:
            edges_t = np.sign(edges) * np.log1p(np.abs(edges))
    else:
        edges_t = edges.copy()

    norm = BoundaryNorm(edges_t, ncolors=n_bands)
    return ColorScheme(
        mode=mode, edges=edges, edges_t=edges_t, cmap=cmap, norm=norm, log=log, cap=cap
    )
