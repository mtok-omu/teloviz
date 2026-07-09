# teloviz

Telomere repeat detection and visualization along genome sequences.

`teloviz` scans a genome FASTA for telomere repeat motifs (default `TTAGGG`,
searched on both strands), summarizes their density toward sequence ends, and
renders a figure so you can see which contigs/chromosomes are capped by
telomeres and where internal telomeric repeats occur.

> **Status:** early skeleton. The CLI parses inputs; detection and plotting are
> not implemented yet.

## Install (dev)

```bash
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate tool_dev
pip install -e .
```

## Usage

```bash
teloviz --input genome.fasta --motif TTAGGG --window 1000 --output teloviz_out
```

| Option | Default | Meaning |
|--------|---------|---------|
| `-i, --input` | (required) | Genome FASTA (assembly or contigs) |
| `-m, --motif` | `TTAGGG` | Telomere repeat motif (both strands) |
| `-w, --window` | `1000` | Window size (bp) for density binning |
| `-o, --output` | `teloviz_out` | Output directory for figures + summary |

## Development

```bash
pytest
```
