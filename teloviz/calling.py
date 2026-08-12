"""Telomere-cap calling for teloviz (spec section 10).

For each chromosome, sum the telomeric-repeat count (forward + reverse, after the
same --motif/--min-count filtering the plot uses) over the windows within
``dist_bp`` of each end. If that end-region total is >= ``threshold``, the end is
called telomere-capped. The call is a property of the assembly, not of the color
scheme, so ``sum`` and ``orientation`` plots share one set of calls.

The window position is its *end* coordinate (tidk convention), so a window is a
5' end window when ``window <= dist_bp`` and a 3' end window when
``window > length - dist_bp``. tidk windows are typically 10 kb, so ``dist_bp``
is naturally a small multiple of that.
"""

from __future__ import annotations

from dataclasses import dataclass

from .prepare import Prepared

# Below this multiple of the call distance the 5' and 3' end regions overlap.
SHORT_FACTOR = 2


def is_short(length: int, dist_bp: int) -> bool:
    """True if ``length`` is too short for the two end regions to be disjoint.

    The 5' region is ``window <= dist_bp`` and the 3' region is
    ``window > length - dist_bp``, so below ``SHORT_FACTOR * dist_bp`` they
    overlap and the same windows are counted for both ends — a sequence with
    telomeric repeats at only one physical end can then be called capped at
    both. Single source of truth for the plot and the HTML report.
    """
    return dist_bp > 0 and length < SHORT_FACTOR * dist_bp


@dataclass
class Call:
    """Telomere-cap verdict for one chromosome, both ends."""

    id: str
    length: int
    five_count: int      # forward+reverse summed over the 5' end region
    three_count: int     # forward+reverse summed over the 3' end region
    five: bool           # 5' end called capped (five_count >= threshold)
    three: bool          # 3' end called capped
    short: bool = False  # end regions overlap at this length (call may be spurious)

    @property
    def both(self) -> bool:
        return self.five and self.three

    @property
    def suffix(self) -> str:
        """Name suffix: ``*`` for a both-ends call, ``(*)`` when it may be spurious.

        Parenthesizing (rather than only recoloring the plot) keeps the warning
        readable in grayscale print and for color-vision-deficient readers, and
        keeps the figure and the HTML report labelled the same way.
        """
        if not self.both:
            return ""
        return " (*)" if self.short else " *"

    @property
    def status(self) -> str:
        if self.five and self.three:
            return "both"
        if self.five:
            return "5' only"
        if self.three:
            return "3' only"
        return "none"


def call_telomeres(prepared: Prepared, *, dist_bp: int, threshold: int) -> list[Call]:
    """Return a per-chromosome telomere-cap call, in the plotted (natural) order."""
    table = prepared.table
    calls: list[Call] = []
    for cid in prepared.order:
        length = prepared.lengths[cid]
        sub = table[table["id"] == cid]
        w = sub["window"].to_numpy()
        tot = sub["forward"].to_numpy() + sub["reverse"].to_numpy()
        five_count = int(tot[w <= dist_bp].sum())
        three_count = int(tot[w > length - dist_bp].sum())
        calls.append(Call(
            id=cid, length=length,
            five_count=five_count, three_count=three_count,
            five=five_count >= threshold, three=three_count >= threshold,
            short=is_short(length, dist_bp),
        ))
    return calls
