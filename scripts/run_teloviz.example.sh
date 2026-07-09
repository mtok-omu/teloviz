#!/usr/bin/env bash
# Example invocation. For large genomes, submit via Slurm rather than running on
# the login node.
set -euo pipefail

source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate tool_dev

teloviz \
  --input demo/example.fasta \
  --motif TTAGGG \
  --window 1000 \
  --output teloviz_out
