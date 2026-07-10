#!/usr/bin/env bash
# Example invocation. tidk must already have been run to produce the windows TSV.
# teloviz itself is light and headless (no display, no external binary).
set -euo pipefail

# Activate whatever environment teloviz is installed in, e.g.:
#   conda activate teloviz   # or: source .venv/bin/activate

teloviz \
  sample_telomeric_repeat_windows.tsv \
  --fai sample.fa.fai \
  --mode both \
  --bin 100 --cap 500 \
  -o sample \
  --format pdf,png
