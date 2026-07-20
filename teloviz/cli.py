"""Command-line interface for teloviz (spec section 8).

Orchestrates the pipeline: load tidk TSV (+ optional .fai) -> aggregate/prepare
-> per-mode color scheme -> matplotlib ideogram -> PDF/PNG/SVG export.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__

VALID_FORMATS = ("pdf", "png", "svg")


def parse_formats(value: str) -> list[str]:
    """Parse a comma-separated --format value into a validated ordered list."""
    fmts = [f.strip().lower() for f in value.split(",") if f.strip()]
    bad = [f for f in fmts if f not in VALID_FORMATS]
    if bad:
        raise argparse.ArgumentTypeError(
            f"unknown format(s): {', '.join(bad)} (choose from {', '.join(VALID_FORMATS)})"
        )
    if not fmts:
        raise argparse.ArgumentTypeError("no output format given")
    return fmts


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="teloviz",
        description="Ideogram-style visualization of tidk telomere-repeat windows.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"teloviz {__version__}")
    p.add_argument(
        "input", type=Path, metavar="INPUT",
        help="tidk search '*_telomeric_repeat_windows.tsv'.",
    )
    p.add_argument(
        "-o", "--out-prefix", default=None,
        help="Output filename prefix (default: derived from INPUT name).",
    )
    p.add_argument(
        "--fai", type=Path, default=None,
        help="samtools faidx '.fai' for exact chromosome lengths.",
    )
    p.add_argument(
        "--mode", choices=("sum", "orientation", "both"), default="sum",
        help="sum (forward+reverse) / orientation (forward-reverse, diverging) / both.",
    )
    p.add_argument("--bin", type=int, default=100, help="Color bin width.")
    p.add_argument("--cap", type=float, default=500, help="Color cap; values above clamp.")
    p.add_argument("--log", action="store_true", help="Log-scale color mapping.")
    p.add_argument(
        "--motif", default=None,
        help="Limit to one motif (default: aggregate all motifs).",
    )
    p.add_argument("--cmap-sum", default="Reds", help="Colormap for sum mode.")
    p.add_argument("--cmap-div", default="RdBu_r", help="Diverging colormap for orientation mode.")
    p.add_argument(
        "--style", choices=("dot", "rect"), default="dot",
        help="Mark style on the chromosome bar: dot (fixed-size round marker at "
             "each window, always visible) / rect (honest length-proportional).",
    )
    p.add_argument(
        "--dot-size", type=float, default=40.0,
        help="Round marker area in points^2 (dot style only; larger = bigger dots).",
    )
    p.add_argument(
        "--call-dist", type=float, default=30.0, metavar="KB",
        help="Telomere call: look within this many kb of each chromosome end "
             "(tidk windows are ~10 kb, so this is a small multiple of that).",
    )
    p.add_argument(
        "--call-min", type=int, default=50,
        help="Telomere call: forward+reverse summed over the end region must reach "
             "this to call the end telomere-capped.",
    )
    p.add_argument(
        "--no-call", action="store_true",
        help="Disable telomere calling (no end markers on the plot, no HTML report).",
    )
    p.add_argument(
        "--rDNA", dest="rdna", type=Path, default=None, metavar="BED",
        help="Optional feature BED (rDNA/NOR arrays). Off by default = telomere "
             "only. Drawn as neutral lanes below each bar (one lane per BED 'type' "
             "column) and used to annotate uncapped ends in the report. Expected "
             "pre-merged (bedtools merge); teloviz only visualizes it.",
    )
    p.add_argument(
        "--proximity", type=float, default=500.0, metavar="KB",
        help="Report-only: an uncapped end is annotated with a feature within "
             "this many kb of it (annotation distance, not a QC/color threshold). "
             "Default 500: a real NOR can sit well inside the end (often >100 kb).",
    )
    p.add_argument(
        "--min-len", type=int, default=0,
        help="Drop sequences shorter than this (0 = keep all).",
    )
    p.add_argument(
        "--min-count", type=int, default=0,
        help="Noise floor: windows with forward+reverse below this are treated "
             "as background/white (0 = keep all). Suppresses the pervasive "
             "random-match background so real telomere arrays stand out.",
    )
    p.add_argument("--font-size", type=float, default=10.0, metavar="PT",
                   help="Base text size in points; all labels/title/ticks scale "
                        "from it (default: 10).")
    p.add_argument("--width", type=float, default=None, metavar="MM",
                   help="Figure width in millimeters (default: auto = 254 mm / 10 in).")
    p.add_argument("--height", type=float, default=None, metavar="MM",
                   help="Figure height in millimeters (default: auto, scales with the "
                        "number of chromosomes). Set both --width and --height to "
                        "fix the exact aspect ratio (whitespace trimming is turned "
                        "off so the ratio is preserved).")
    p.add_argument("--dpi", type=int, default=200, help="Raster (PNG) resolution.")
    p.add_argument(
        "--format", type=parse_formats, default="pdf", dest="formats",
        help="Comma-separated output formats: pdf / png / svg.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    # The string default bypasses type=parse_formats (only CLI-supplied values are
    # converted), so normalize it here.
    if isinstance(args.formats, str):
        args.formats = parse_formats(args.formats)

    if not args.input.exists():
        print(f"teloviz: input not found: {args.input}", file=sys.stderr)
        return 2
    if args.fai is not None and not args.fai.exists():
        print(f"teloviz: fai not found: {args.fai}", file=sys.stderr)
        return 2
    if args.rdna is not None and not args.rdna.exists():
        print(f"teloviz: feature BED not found: {args.rdna}", file=sys.stderr)
        return 2

    out_prefix = args.out_prefix or args.input.name.split("_telomeric_repeat_windows")[0]

    # Import the pipeline lazily so `--version`/`--help` don't pay the matplotlib
    # import cost.
    from .calling import call_telomeres
    from .color import build_scheme
    from .features import load_features, normalize_features
    from .io_windows import TelovizInputError, load_fai, load_windows
    from .prepare import prepare
    from .render import render, save
    from .report import write_report
    from ._mpl import plt

    try:
        df = load_windows(args.input)
    except TelovizInputError as e:
        print(f"teloviz: {e}", file=sys.stderr)
        return 2

    fai = load_fai(args.fai) if args.fai is not None else None
    prepared = prepare(
        df, fai=fai, motif=args.motif, min_len=args.min_len, min_count=args.min_count
    )

    if not prepared.order:
        print("teloviz: no sequences to plot (check --motif / --min-len).", file=sys.stderr)
        return 1

    # Optional feature track (off by default). Normalized once against the plotted
    # chromosome lengths; independent of any color cutoff.
    features = None
    if args.rdna is not None:
        features = normalize_features(load_features(args.rdna), prepared.lengths)
    proximity_bp = int(args.proximity * 1000)

    # Telomere calling is one verdict per assembly (mode-independent). Shown on
    # every plot and, unless disabled, exported as a standalone HTML report.
    dist_bp = int(args.call_dist * 1000)
    if args.no_call:
        calls, call_note = None, None
    else:
        calls = call_telomeres(prepared, dist_bp=dist_bp, threshold=args.call_min)
        call_note = (f"▸◂ telomere call: forward+reverse ≥ {args.call_min} "
                     f"within {args.call_dist:g} kb of an end   ·   * = both ends")

    modes = ["sum", "orientation"] if args.mode == "both" else [args.mode]
    written: list[str] = []
    for mode in modes:
        scheme = build_scheme(
            mode, bin_w=args.bin, cap=args.cap, log=args.log,
            cmap_sum=args.cmap_sum, cmap_div=args.cmap_div,
        )
        # Figure sizes are given in mm on the CLI; matplotlib is inch-native.
        width_in = args.width / 25.4 if args.width is not None else None
        height_in = args.height / 25.4 if args.height is not None else None
        fig = render(prepared, scheme, style=args.style, dot_size=args.dot_size,
                     calls=calls, call_note=call_note, font_size=args.font_size,
                     width=width_in, height=height_in, features=features)
        label = mode
        # When the user pins both dimensions, keep the exact size (no tight-bbox
        # trim) so the requested aspect ratio is honored.
        exact_size = args.width is not None and args.height is not None
        paths = save(fig, out_prefix, label, args.formats, args.dpi,
                     tight=not exact_size)
        plt.close(fig)
        written.extend(str(p) for p in paths)

    if calls is not None:
        report_path = write_report(
            f"{out_prefix}.telomere_report.html", calls,
            meta={
                "command": "teloviz " + " ".join(argv if argv is not None else sys.argv[1:]),
                "input": str(args.input),
                "mode": args.mode,
                "dist_kb": f"{args.call_dist:g}",
                "call_min": args.call_min,
                "min_count": args.min_count,
                "motif": args.motif,
                "window_size": prepared.window_size,
                "features": str(args.rdna) if args.rdna is not None else "",
            },
            features=features,
            proximity_bp=proximity_bp,
        )
        written.append(str(report_path))

    print(
        f"teloviz {__version__}: {len(prepared.order)} sequences, "
        f"window_size={prepared.window_size}. Wrote:"
    )
    for p in written:
        print(f"  {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
