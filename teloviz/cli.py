"""Command-line interface for teloviz (spec section 8).

Skeleton scope: parse + validate arguments (the full CLI surface from the spec).
The data pipeline (load tidk TSV -> aggregate -> bin to colors -> matplotlib
ideogram -> PDF/PNG/SVG export) is wired in later steps.
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

    # TODO: pipeline — load(args.input) -> aggregate(args.motif) -> bin(args.bin,
    # args.cap, args.log) -> matplotlib ideogram per mode -> savefig(args.formats).
    print(
        f"teloviz {__version__}: parsed OK "
        f"(input={args.input}, mode={args.mode}, bin={args.bin}, cap={args.cap}, "
        f"formats={','.join(args.formats)}, out_prefix={out_prefix}). "
        f"Pipeline not yet implemented.",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
