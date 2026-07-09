# teloviz

Ideogram-style visualization of telomere repeats for assembly QC.

`teloviz` reads the per-window telomere-repeat table from
[`tidk search`](https://github.com/tolkit/telomeric-identifier)
(`*_telomeric_repeat_windows.tsv`), draws each chromosome as a single
length-proportional horizontal bar, and color-codes every window by repeat
amount. Unlike `tidk plot` (per-window line charts), it lays all chromosomes
side by side so end-to-end capping, strand orientation, and internal
(interstitial) repeat clusters — a misjoin signal — are visible at a glance.

Two complementary modes:
- **sum** (default) — `forward + reverse`, white→red. Catches balanced internal
  clusters the orientation mode can miss.
- **orientation** — `forward − reverse` (net), diverging blue←white→red. Shows
  strand structure: normal ends flip red at one tip / blue at the other.

Rendered with **matplotlib** to vector **PDF** (default; also PNG/SVG). Fully
headless — no browser or external binary needed.

> **Status:** working v0.1 — full pipeline (load → prepare → color → render) for
> `sum` / `orientation` / `both`, PDF/PNG/SVG output. Try it on the bundled demo:
> `python demo/make_demo.py && teloviz demo/demo_telomeric_repeat_windows.tsv --fai demo/demo.fa.fai --mode both -o demo/sample`.
> Spec: `tool_dev.md` (internal).

## Install (dev)

```bash
source $HOME/miniforge3/etc/profile.d/conda.sh
conda activate tool_dev
pip install -e .
```

## Usage

```bash
# default: sum mode, white→red, bin 100 / cap 500, PDF
teloviz sample_telomeric_repeat_windows.tsv --fai sample.fa.fai -o sample

# orientation mode, tuned bin/cap, log scale
teloviz sample_..._windows.tsv --mode orientation --bin 50 --cap 1000 --log -o sample_orient

# both modes at once
teloviz sample_..._windows.tsv --mode both -o sample
```

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
| `--min-len` | `0` | Drop shorter sequences (0 = keep all) |
| `--width` / `--height` | auto | Figure size (px) |
| `--dpi` | `200` | Raster (PNG) resolution |
| `--format` | `pdf` | `pdf` / `png` / `svg`, comma-separated |

## Development

```bash
pytest
```
