# teloviz

Ideogram-style visualization of telomere repeats for genome-assembly QC.

`teloviz` takes the per-window telomere-repeat table produced by
[`tidk search`](https://github.com/tolkit/telomeric-identifier) and draws every
chromosome as a length-proportional bar, marking where telomeric repeats occur.
Unlike `tidk plot` (one line chart per sequence), it lays all chromosomes side
by side, so you can see at a glance:

- **end-to-end capping** — is each chromosome telomere-capped at both ends (T2T)?
- **strand orientation** — do the ends show the normal forward/reverse flip?
- **interstitial clusters** — telomeric repeats *inside* a chromosome, a misjoin signal.

It also **calls** which ends are telomere-capped and writes a summary HTML report.

Optionally, an **annotation track** (`--rDNA`, off by default) draws feature
intervals — rDNA/NOR arrays — as a color-coded lane *below* each bar and
annotates un-capped ends in the report (e.g. *cap missing — rDNA at this end*),
so you can tell a genuine missing telomere from an assembly that simply
**stopped at the NOR**. The bar colors are never touched.

Output is vector **PDF** (default; also PNG/SVG), rendered with matplotlib and
fully headless — no browser or external binary needed, so it runs on an HPC
login/compute node as-is.

---

## Requirements

- **Python ≥ 3.10**
- Python packages `pandas`, `matplotlib`, `natsort` (installed automatically below)
- **[tidk](https://github.com/tolkit/telomeric-identifier)** — *upstream, run separately.*
  teloviz only visualizes tidk's output; it does not detect repeats itself.
- Optional: `samtools` (to make the `.fai` index for exact chromosome lengths)

## Installation

```bash
git clone git@github.com:mtok-omu/teloviz.git
cd teloviz

python -m venv .venv && source .venv/bin/activate   # or use a conda env
pip install .
```

This installs the `teloviz` command and its dependencies. Check it:

```bash
teloviz --version
teloviz --help      # full list of options
```

## Quick start

Verify your install on the bundled synthetic demo (no tidk needed):

```bash
python demo/make_demo.py
teloviz demo/demo_telomeric_repeat_windows.tsv --fai demo/demo.fa.fai --mode both -o demo/out
# -> demo/out.sum.pdf, demo/out.orientation.pdf, demo/out.telomere_report.html
```

## Usage on your own assembly

teloviz is step 2 of a two-step workflow — run tidk first, then teloviz:

```bash
# 1. Detect telomeric repeats per window with tidk (see tidk's docs for flags).
#    This produces  <name>_telomeric_repeat_windows.tsv
tidk search --string TTAGGG --output myasm --dir tidk_out genome.fasta

# 2. (optional) exact chromosome lengths for correct bar lengths / end clamping
samtools faidx genome.fasta          # -> genome.fasta.fai

# 3. Visualize + call telomeres with teloviz
teloviz tidk_out/myasm_telomeric_repeat_windows.tsv \
    --fai genome.fasta.fai --mode both -o myasm
```

More examples:

```bash
# suppress random-match background (recommended on real data), bigger dots
teloviz windows.tsv --fai genome.fa.fai --min-count 50 --dot-size 60 -o myasm

# orientation mode with a tuned scale and log color mapping
teloviz windows.tsv --mode orientation --bin 50 --cap 1000 --log -o myasm

# honest length-proportional rectangles instead of fixed-size dots
teloviz windows.tsv --style rect -o myasm

# fix the exact aspect ratio (inches); e.g. a wide 16x6 figure
teloviz windows.tsv --width 16 --height 6 -o myasm

# add an annotation track (rDNA/NOR arrays) below the bars
teloviz windows.tsv --fai genome.fa.fai --rDNA features.bed -o myasm
```

### Annotation track (`--rDNA`)

`--rDNA features.bed` overlays genomic features and is **off by default** (no BED
→ telomere only). The BED is expected to be pre-merged (`bedtools merge`);
teloviz only visualizes it — it does not cluster or count features.

Tab-separated, no header, `#` comments; **columns**: `chrom  start  end  [name]
[type]` (BED coordinates: 0-based start, exclusive end). The 5th column `type`
sets the **color** — every distinct value gets its own color on one shared lane
below the bar (a legend maps color → type). Adding a new `type` just adds a
color, no code change (future `centromere` / `gap`):

```
chr7    12300000   12520000   45S_n28   rdna_45S
chr4    8000000    8330000    5S_n189   rdna_5S
```

In the report, each **un-capped** end gets a note naming the nearest feature
within `--proximity` kb (default 500), or stating *no feature nearby* — because
"no rDNA and still no telomere" is itself evidence for a real gap. (500 kb by
default because a real NOR can sit well inside the end, often >100 kb in.)

**Building the BED from a barrnap GFF.** teloviz reads BED only (it does not
parse GFF). If your features come from `barrnap`/`pybarrnap`, convert once with
the bundled helper — a preprocessing step, run before teloviz:

```bash
scripts/rrna_gff_to_bed.sh rRNA.gff features.bed          # 45S NOR arrays
scripts/rrna_gff_to_bed.sh --with-5s rRNA.gff feats.bed   # also add the 5S array
scripts/rrna_gff_to_bed.sh --dry-run rRNA.gff             # show what the cutoffs do
```

For `rdna_45S` it merges 18S/5.8S/28S hits and keeps an array only where **all
three components are present** — a complete 45S unit is 18S-ITS1-5.8S-ITS2-28S,
so this is a threshold-free binary test that separates a real NOR from a stray
18S/28S fragment. `--with-5s` adds the 5S array (a single gene, so it uses a
copy-count cutoff). Needs `bedtools`. The chromosome names it emits (from the
GFF) must match your tidk TSV / `.fai`.

## Reading the plot

- Each chromosome is one gray length bar; every window with telomeric repeats is
  a colored dot at its true position (`--style dot`, default — a fixed on-screen
  size so tiny telomere windows stay visible; `--style rect` draws true-width
  rectangles instead).
- **sum** mode (default): `forward + reverse`, white→red. Catches balanced
  internal clusters that orientation mode can miss.
- **orientation** mode: `forward − reverse` (net), diverging blue←white→red.
  Normal ends flip red at one tip and blue at the other. `--mode both` writes both.
- **Telomere call:** an end is called capped when `forward+reverse` within
  `--call-dist` kb of it reaches `--call-min` (defaults 30 kb / 50). Capped ends
  get an inward triangle (5′ ▸ / 3′ ◂); a chromosome capped at *both* ends gets a
  `*`. Disable with `--no-call`.

## Output files

For `-o myasm` (and `--format pdf,png` to add PNGs):

| File | What it is |
|------|-----------|
| `myasm.sum.pdf` / `myasm.orientation.pdf` | the ideogram(s) for the chosen mode(s) |
| `myasm.telomere_report.html` | per-chromosome end counts, calls, and the exact settings/command used |

## Options

`teloviz --help` is the authoritative, always-current list. Summary:

| Option | Default | Meaning |
|--------|---------|---------|
| `INPUT` | (required) | tidk `*_telomeric_repeat_windows.tsv` |
| `-o, --out-prefix` | input name | Output filename prefix |
| `--fai` | none | `.fai` for exact chromosome lengths |
| `--mode` | `sum` | `sum` / `orientation` / `both` |
| `--bin` | `100` | Color bin width |
| `--cap` | `500` | Color cap (values above clamp) |
| `--log` | off | Log-scale color mapping |
| `--motif` | none | Limit to one motif (default: all) |
| `--cmap-sum` / `--cmap-div` | `Reds` / `RdBu_r` | Colormaps |
| `--style` | `dot` | Mark style: `dot` (fixed-size marker, always visible) / `rect` (length-proportional) |
| `--dot-size` | `40` | Round marker area in points² (dot style only) |
| `--call-dist` | `30` | Telomere call: look within this many kb of each end |
| `--call-min` | `50` | Telomere call: forward+reverse in that end region ≥ this → capped |
| `--no-call` | off | Disable telomere calling (no end markers, no HTML report) |
| `--rDNA` | none | Optional feature BED (rDNA/NOR); lanes below bars + report annotation |
| `--proximity` | `500` | Report-only: annotate an un-capped end with a feature within this many kb |
| `--min-len` | `0` | Drop sequences shorter than this (0 = keep all) |
| `--min-count` | `0` | Noise floor: windows with forward+reverse below this → white (both modes) |
| `--font-size` | `10` | Base text size in points; all labels/title/ticks scale from it |
| `--width` / `--height` | auto | Figure size in **inches**. Set both to fix the exact aspect ratio (whitespace trim is disabled so the ratio is preserved) |
| `--dpi` | `200` | Raster (PNG) resolution (output px = inches × dpi) |
| `--format` | `pdf` | `pdf` / `png` / `svg`, comma-separated |

## Development

```bash
pip install -e .[test]   # editable install + pytest
pytest
```
