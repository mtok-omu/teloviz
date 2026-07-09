"""Input readers for teloviz (spec section 3).

- ``load_windows``: the required tidk ``*_telomeric_repeat_windows.tsv``.
- ``load_fai``: an optional ``samtools faidx`` ``.fai`` for exact lengths.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

# tidk column -> internal column
_COLMAP = {
    "id": "id",
    "window": "window",
    "forward_repeat_number": "forward",
    "reverse_repeat_number": "reverse",
    "telomeric_repeat": "motif",
}


class TelovizInputError(ValueError):
    """Raised when an input file is missing required columns or malformed."""


def load_windows(path: Path) -> pd.DataFrame:
    """Load a tidk windows TSV into columns id/window/forward/reverse/motif."""
    df = pd.read_csv(path, sep="\t")
    missing = [c for c in _COLMAP if c not in df.columns]
    if missing:
        raise TelovizInputError(
            f"{path}: missing required column(s): {', '.join(missing)} "
            f"(found: {', '.join(map(str, df.columns))})"
        )
    df = df[list(_COLMAP)].rename(columns=_COLMAP)
    for col in ("window", "forward", "reverse"):
        df[col] = pd.to_numeric(df[col], errors="raise").astype("int64")
    df["id"] = df["id"].astype(str)
    df["motif"] = df["motif"].astype(str)
    return df


def load_fai(path: Path) -> dict[str, int]:
    """Load sequence lengths from a ``.fai`` (columns: name, length, ...)."""
    fai = pd.read_csv(path, sep="\t", header=None, usecols=[0, 1], names=["id", "length"])
    return {str(k): int(v) for k, v in zip(fai["id"], fai["length"])}
