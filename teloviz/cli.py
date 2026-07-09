"""Command-line interface for teloviz.

Skeleton scope: parse arguments and validate inputs. The detection and plotting
stages (scan telomere repeats -> summarize per-sequence ends -> render figure)
are wired in later once the spec is settled.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__

# Default telomere repeat motif. Vertebrates use TTAGGG; many plants use
# TTTAGGG. Overridable with --motif.
DEFAULT_MOTIF = "TTAGGG"


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="teloviz",
        description="Detect and visualize telomere repeats along genome sequences.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--version", action="version", version=f"teloviz {__version__}")
    p.add_argument(
        "-i", "--input", type=Path, required=True,
        help="Input genome FASTA (assembly or contigs).",
    )
    p.add_argument(
        "-m", "--motif", default=DEFAULT_MOTIF,
        help="Telomere repeat motif to search for (both strands).",
    )
    p.add_argument(
        "-w", "--window", type=int, default=1000,
        help="Window size (bp) for repeat-density binning.",
    )
    p.add_argument(
        "-o", "--output", type=Path, default=Path("teloviz_out"),
        help="Output directory for the figure(s) and summary table.",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if not args.input.exists():
        print(f"teloviz: input not found: {args.input}", file=sys.stderr)
        return 2

    # TODO: scan(args.input, args.motif, args.window) -> render(args.output)
    print(
        f"teloviz {__version__}: parsed OK "
        f"(input={args.input}, motif={args.motif}, window={args.window}, "
        f"output={args.output}). Detection/plotting not yet implemented.",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
