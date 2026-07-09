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
        "--min-len", type=int, default=0,
        help="Drop sequences shorter than this (0 = keep all).",
    )
    p.add_argument(
        "--min-count", type=int, default=0,
        help="Noise floor: windows with forward+reverse below this are treated "
             "as background/white (0 = keep all). Suppresses the pervasive "
             "random-match background so real telomere arrays stand out.",
    )
    p.add_argument("--width", type=int, default=None, help="Figure width in px (auto).")
    p.add_argument("--height", type=int, default=None, help="Figure height in px (auto).")
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

    out_prefix = args.out_prefix or args.input.name.split("_telomeric_repeat_windows")[0]

    # Import the pipeline lazily so `--version`/`--help` don't pay the matplotlib
    # import cost.
    from .color import build_scheme
    from .io_windows import TelovizInputError, load_fai, load_windows
    from .prepare import prepare
    from .render import render, save
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

    modes = ["sum", "orientation"] if args.mode == "both" else [args.mode]
    written: list[str] = []
    for mode in modes:
        scheme = build_scheme(
            mode, bin_w=args.bin, cap=args.cap, log=args.log,
            cmap_sum=args.cmap_sum, cmap_div=args.cmap_div,
        )
        fig = render(prepared, scheme, width=args.width, height=args.height)
        paths = save(fig, out_prefix, mode, args.formats, args.dpi)
        plt.close(fig)
        written.extend(str(p) for p in paths)

    print(
        f"teloviz {__version__}: {len(prepared.order)} sequences, "
        f"window_size={prepared.window_size}. Wrote:"
    )
    for p in written:
        print(f"  {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
