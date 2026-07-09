"""Generate a small synthetic tidk-style windows TSV + .fai for demos/tests.

Deterministic (no RNG). Models a few chromosomes with:
  - a normal telomere: forward-dominant at the 5' end, reverse-dominant at 3'
  - one chromosome with an internal (interstitial) repeat cluster = misjoin hint
Run:  python demo/make_demo.py
"""

from __future__ import annotations

from pathlib import Path

WINDOW = 10000
HEADER = "id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\n"

# id -> (n_windows, internal_cluster_window_index or None)
CHROMS = {
    "chr1": (60, None),
    "chr2": (45, 22),   # internal cluster ~ midpoint (misjoin suspect)
    "chr10": (80, None),
}


def repeats(k: int, n: int, cluster: int | None) -> tuple[int, int]:
    """(forward, reverse) counts for window k (1-based) of an n-window chrom."""
    fwd = rev = 0
    if k <= 2:                    # 5' telomere: forward-dominant
        fwd, rev = 340 - 120 * (k - 1), 3
    if k >= n - 1:                # 3' telomere: reverse-dominant
        rev, fwd = 340 - 120 * (n - k), 3
    if cluster is not None and abs(k - cluster) <= 1:  # balanced internal cluster
        fwd, rev = 180, 170
    return fwd, rev


def main() -> None:
    here = Path(__file__).resolve().parent
    lines = [HEADER]
    fai = []
    for cid, (n, cluster) in CHROMS.items():
        for k in range(1, n + 1):
            fwd, rev = repeats(k, n, cluster)
            lines.append(f"{cid}\t{k * WINDOW}\t{fwd}\t{rev}\tTTAGGG\n")
        fai.append(f"{cid}\t{n * WINDOW}\t0\t{WINDOW}\t{WINDOW + 1}\n")

    (here / "demo_telomeric_repeat_windows.tsv").write_text("".join(lines))
    (here / "demo.fa.fai").write_text("".join(fai))
    print("wrote demo_telomeric_repeat_windows.tsv and demo.fa.fai")


if __name__ == "__main__":
    main()
