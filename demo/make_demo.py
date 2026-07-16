"""Generate a small synthetic tidk-style windows TSV + .fai (+ features BED) for
demos/tests.

Deterministic (no RNG). Models the situations the annotation track exists to
tell apart (annotation-track spec section 6):
  - normal T2T telomere: forward-dominant 5' end, reverse-dominant 3' end
  - one chromosome with an internal (interstitial) repeat cluster = misjoin hint
  - three feature scenarios (all rdna_45S, one shared lane):
    1. chr7  = NOR-stop / "false loss": 3' end un-capped + a single 45S unit
       right at that end  -> the telomere is probably there but the assembly
       stopped in the array. THE REASON THIS TRACK EXISTS.
    2. chr8  = penetrated: 3' end capped AND a multi-unit 45S array sitting well
       inside it, with sequence beyond -> rDNA present yet the cap is normal.
       The contrast that keeps scenario 1 honest.
    3. chr5  = 5' end un-capped with NO feature nearby -> a genuine gap suspect.
Run:  python demo/make_demo.py
"""

from __future__ import annotations

from pathlib import Path

WINDOW = 10000
HEADER = "id\twindow\tforward_repeat_number\treverse_repeat_number\ttelomeric_repeat\n"

# id -> (n_windows, internal_cluster_window_index or None, ends_with_telomere)
CHROMS = {
    "chr1": (60, None, ("5", "3")),
    "chr2": (45, 22, ("5", "3")),    # internal cluster ~ midpoint (misjoin suspect)
    "chr5": (50, None, ("3",)),      # 5' un-capped, no feature nearby
    "chr7": (70, None, ("5",)),      # 3' un-capped, single 45S unit at that end
    "chr8": (90, None, ("5", "3")),  # capped both ends, 45S array inside the 3' end
    "chr10": (80, None, ("5", "3")),
}

# Feature BED (annotation track). 45S only -> a single lane (spec section 2.3).
FEATURES = [
    # chr7 (length 700000): a single 45S unit ~10 kb short of the un-capped 3' end.
    ("chr7", 680000, 690000, "45S_n1", "rdna_45S"),
    # chr8 (length 900000): a 45S array at 650-720 kb, then ~180 kb of sequence to
    # the capped 3' end -> "penetrated".
    ("chr8", 650000, 720000, "45S_n7", "rdna_45S"),
]


def repeats(k: int, n: int, cluster: int | None, ends: tuple[str, ...]) -> tuple[int, int]:
    """(forward, reverse) counts for window k (1-based) of an n-window chrom."""
    fwd = rev = 0
    if "5" in ends and k <= 2:            # 5' telomere: forward-dominant
        fwd, rev = 340 - 120 * (k - 1), 3
    if "3" in ends and k >= n - 1:        # 3' telomere: reverse-dominant
        rev, fwd = 340 - 120 * (n - k), 3
    if cluster is not None and abs(k - cluster) <= 1:  # balanced internal cluster
        fwd, rev = 180, 170
    return fwd, rev


def main() -> None:
    here = Path(__file__).resolve().parent
    lines = [HEADER]
    fai = []
    for cid, (n, cluster, ends) in CHROMS.items():
        for k in range(1, n + 1):
            fwd, rev = repeats(k, n, cluster, ends)
            lines.append(f"{cid}\t{k * WINDOW}\t{fwd}\t{rev}\tTTAGGG\n")
        fai.append(f"{cid}\t{n * WINDOW}\t0\t{WINDOW}\t{WINDOW + 1}\n")

    tsv = here / "demo_telomeric_repeat_windows.tsv"
    fai_path = here / "demo.fa.fai"
    bed_path = here / "demo_features.bed"
    tsv.write_text("".join(lines))
    fai_path.write_text("".join(fai))
    bed_path.write_text(
        "".join(f"{c}\t{s}\t{e}\t{name}\t{t}\n" for c, s, e, name, t in FEATURES)
    )
    print(f"wrote {tsv}")
    print(f"wrote {fai_path}")
    print(f"wrote {bed_path}")
    print("\nTry:")
    print(f"  teloviz {tsv} --fai {fai_path} --mode both \\")
    print(f"      --rDNA {bed_path} -o teloviz_out/demo")


if __name__ == "__main__":
    main()
