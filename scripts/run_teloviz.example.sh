#!/usr/bin/env bash
# Example invocation. tidk must already have been run to produce the windows TSV.
# teloviz itself is light and headless (no display, no external binary).
set -euo pipefail

source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate tool_dev

teloviz \
  sample_telomeric_repeat_windows.tsv \
  --fai sample.fa.fai \
  --mode both \
  --bin 100 --cap 500 \
  -o sample \
  --format pdf,png
