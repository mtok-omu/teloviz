"""Data preparation for teloviz (spec section 4).

Aggregate motifs, estimate a single genome-wide window size, resolve chromosome
lengths (from .fai or max window), drop short sequences, and natural-sort ids.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from natsort import natsorted

DEFAULT_WINDOW_SIZE = 10000


@dataclass
class Prepared:
    """Everything the renderer needs, after cleaning."""

    order: list[str]              # chromosome ids, natural-sorted (top -> bottom)
    lengths: dict[str, int]       # id -> chromosome length (bp)
    table: pd.DataFrame           # columns: id, window, forward, reverse
    window_size: int              # single genome-wide value


def estimate_window_size(table: pd.DataFrame) -> int:
    """Most common spacing of ``window`` end coords across the whole file.

    tidk uses one ``-w`` for every sequence in a run, so the window size is a
    property of the file, not of a single chromosome. We pool the diffs from all
    chromosomes and take the mode. Sequences with a single window contribute no
    diff and simply inherit this value. Falls back to ``DEFAULT_WINDOW_SIZE``
    only when no chromosome has >= 2 windows.
    """
    diffs: list[int] = []
    for _, g in table.groupby("id"):
        w = np.sort(g["window"].unique())
        if w.size >= 2:
            diffs.extend(np.diff(w).tolist())
    if not diffs:
        return DEFAULT_WINDOW_SIZE
    values, counts = np.unique(np.asarray(diffs), return_counts=True)
    return int(values[int(np.argmax(counts))])


def prepare(
    df: pd.DataFrame,
    *,
    fai: dict[str, int] | None = None,
    motif: str | None = None,
    min_len: int = 0,
    min_count: int = 0,
) -> Prepared:
    """Filter/aggregate the raw windows table into a ``Prepared`` bundle."""
    if motif is not None:
        df = df[df["motif"] == motif]

    # Collapse motifs: sum forward/reverse per (id, window).
    agg = (
        df.groupby(["id", "window"], as_index=False)[["forward", "reverse"]]
        .sum()
        .sort_values(["id", "window"])
        .reset_index(drop=True)
    )

    window_size = estimate_window_size(agg)

    # Noise floor: windows whose total repeat count (forward + reverse) is below
    # min_count are treated as background (zeroed -> white) in BOTH modes, so the
    # ubiquitous low-level random-match background does not wash out real arrays.
    if min_count > 0:
        below = (agg["forward"] + agg["reverse"]) < min_count
        agg.loc[below, ["forward", "reverse"]] = 0

    lengths: dict[str, int] = {}
    for cid, g in agg.groupby("id"):
        if fai is not None and cid in fai:
            lengths[cid] = int(fai[cid])
        else:
            lengths[cid] = int(g["window"].max())

    keep = [cid for cid, L in lengths.items() if L >= min_len]
    agg = agg[agg["id"].isin(keep)].reset_index(drop=True)
    lengths = {cid: lengths[cid] for cid in keep}
    order = [str(c) for c in natsorted(keep)]

    return Prepared(order=order, lengths=lengths, table=agg, window_size=window_size)
