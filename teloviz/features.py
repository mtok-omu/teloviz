"""Optional annotation-track features for teloviz (annotation-track spec).

A second, optional input: a BED of genomic-interval features (rDNA/NOR arrays
today; centromere/gap in future) that are drawn as thin lanes *below* each
chromosome bar and used to *modify* the telomere-cap report (e.g. "cap missing —
rDNA at this end"). tidk detection does the telomere work; features only add
context. The chromosome-bar colors are never touched.

Design notes (from the spec):

- One lane per distinct ``type`` value; lanes are generated dynamically, so a BED
  with a new ``type`` grows a new lane with no code change.
- Features do NOT depend on any color cutoff; parsing/normalization happens once.
- Clustering/statistics are out of scope: the BED is expected to be pre-merged
  (``bedtools merge``); teloviz only visualizes and annotates.
"""

from __future__ import annotations

import sys
from dataclasses import dataclass


@dataclass
class Feature:
    """One genomic-interval annotation (0-based, half-open, like BED)."""

    chrom: str
    start: int  # 0-based, inclusive
    end: int    # 0-based, exclusive
    name: str   # tooltip/label text (BED col 4; default "chrom:start-end")
    type: str   # lane key (BED col 5; default "feature")


@dataclass
class FeatureSet:
    """Normalized features plus the stable lane order (type first-appearance)."""

    features: list[Feature]     # kept features, after normalization
    lane_order: list[str]       # distinct type values, first-appearance order

    def for_chrom(self, cid: str) -> list[Feature]:
        return [f for f in self.features if f.chrom == cid]

    def __bool__(self) -> bool:
        return bool(self.features)


# Hardcoded type -> color. Deliberately avoids the blue-white-red diverging range
# the orientation bars use, so features never read as telomere signal.
TYPE_COLORS: dict[str, str] = {
    "rdna_45S": "#2E7D32",   # green
    "rdna_5S": "#9E9E9E",    # gray
    "centromere": "#6A1B9A",  # purple (future)
    "gap": "#455A64",         # blue-gray (future)
}

# Fallback palette for unknown types (assigned in first-appearance order). Also
# steers clear of pure blue/red so it cannot be confused with the bar heat.
_PALETTE = ["#00897B", "#F9A825", "#8D6E63", "#5E35B1", "#C2185B", "#546E7A"]

# Types whose presence is a prediction rather than a measurement get a
# "(candidate)" label suffix (spec section 4.3). rDNA from pybarrnap is measured,
# so it is not listed here; centromere predictions would be.
CANDIDATE_TYPES: frozenset[str] = frozenset({"centromere"})

FEATURE_ALPHA = 0.6


def color_for_type(ftype: str, lane_order: list[str]) -> str:
    """Color for a lane ``type``: the hardcoded map, else palette by lane index."""
    if ftype in TYPE_COLORS:
        return TYPE_COLORS[ftype]
    unknown = [t for t in lane_order if t not in TYPE_COLORS]
    idx = unknown.index(ftype) if ftype in unknown else 0
    return _PALETTE[idx % len(_PALETTE)]


def lane_label(ftype: str) -> str:
    """Legend/lane label for a type; predicted types get a '(candidate)' suffix."""
    return f"{ftype} (candidate)" if ftype in CANDIDATE_TYPES else ftype


def load_features(path) -> list[Feature]:
    """Parse a BED-like feature file into raw ``Feature`` records.

    Tab-separated, no header, ``#`` lines are comments. At least 3 columns
    (chrom, start, end); column 4 (name) defaults to ``chrom:start-end`` and
    column 5 (type) defaults to ``feature``. Coordinates are BED (0-based start,
    exclusive end) and passed through unchanged.
    """
    feats: list[Feature] = []
    with open(path, encoding="utf-8") as fh:
        for lineno, raw in enumerate(fh, 1):
            line = raw.rstrip("\n")
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            cols = line.split("\t")
            if len(cols) < 3:
                raise ValueError(
                    f"{path}:{lineno}: need >= 3 tab-separated columns "
                    f"(chrom, start, end); got {len(cols)}"
                )
            chrom = cols[0]
            try:
                start, end = int(cols[1]), int(cols[2])
            except ValueError as exc:
                raise ValueError(
                    f"{path}:{lineno}: start/end must be integers ({exc})"
                ) from None
            name = cols[3] if len(cols) >= 4 and cols[3] != "" else f"{chrom}:{start}-{end}"
            ftype = cols[4] if len(cols) >= 5 and cols[4] != "" else "feature"
            feats.append(Feature(chrom=chrom, start=start, end=end, name=name, type=ftype))
    return feats


def normalize_features(
    raw: list[Feature], lengths: dict[str, int], *, quiet: bool = False
) -> FeatureSet:
    """Filter/clamp features against the plotted chromosomes.

    - Drop features whose ``chrom`` is not among ``lengths`` (report only a count
      to stderr, not one line each).
    - Drop ``end <= start``.
    - Clamp ``end`` to the chromosome length; drop if that leaves nothing.
    - Lane order is the first-appearance order of ``type`` among kept features.
    """
    kept: list[Feature] = []
    lane_order: list[str] = []
    n_unknown_chrom = 0
    n_bad_interval = 0
    for f in raw:
        if f.chrom not in lengths:
            n_unknown_chrom += 1
            continue
        if f.end <= f.start:
            n_bad_interval += 1
            continue
        length = lengths[f.chrom]
        start = max(0, f.start)
        end = min(f.end, length)
        if end <= start:
            n_bad_interval += 1
            continue
        kept.append(Feature(chrom=f.chrom, start=start, end=end, name=f.name, type=f.type))
        if f.type not in lane_order:
            lane_order.append(f.type)

    if not quiet:
        if n_unknown_chrom:
            print(
                f"teloviz: warning: dropped {n_unknown_chrom} feature(s) on "
                f"sequences absent from the tidk/fai data.",
                file=sys.stderr,
            )
        if n_bad_interval:
            print(
                f"teloviz: warning: dropped {n_bad_interval} feature(s) with an "
                f"empty interval (end <= start, or fully past the chromosome end).",
                file=sys.stderr,
            )
    return FeatureSet(features=kept, lane_order=lane_order)
