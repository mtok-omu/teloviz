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


@dataclass
class Call:
    """Telomere-cap verdict for one chromosome, both ends."""

    id: str
    length: int
    five_count: int      # forward+reverse summed over the 5' end region
    three_count: int     # forward+reverse summed over the 3' end region
    five: bool           # 5' end called capped (five_count >= threshold)
    three: bool          # 3' end called capped

    @property
    def both(self) -> bool:
        return self.five and self.three

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
        ))
    return calls
