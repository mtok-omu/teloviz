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

**Telomere calling:** each end is called telomere-capped when `forward+reverse`
within `--call-dist` kb of it reaches `--call-min` (defaults 30 kb / 50). Capped
ends get an inward triangle (5′ ▸ / 3′ ◂) on the plot, both-ends-capped
chromosomes get a `*`, and a standalone `<prefix>.telomere_report.html` lists the
per-chromosome counts, calls, and the exact settings/command (disable: `--no-call`).

> **Status:** working v0.1 — full pipeline (load → prepare → color → render) for
> `sum` / `orientation` / `both`, PDF/PNG/SVG output. Try it on the bundled demo:
> `python demo/make_demo.py && teloviz demo/demo_telomeric_repeat_windows.tsv --fai demo/demo.fa.fai --mode both -o demo/out`.

## Install (dev)

```bash
python -m venv .venv && source .venv/bin/activate   # or activate a conda env
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

# larger dots on a thinner bar
teloviz sample_..._windows.tsv --fai sample.fa.fai --min-count 50 --dot-size 60 -o sample
```

Each non-white window is drawn as a fixed-size dot at its genomic centre on a
length-proportional bar (`--style dot`, default), so telomere arrays — which sit
in a handful of tiny windows — stay visible instead of vanishing to slivers.
Use `--style rect` for honest length-proportional rectangles (zoomable vector).

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
| `--min-count` | `0` | Noise floor: windows with forward+reverse below this → white (both modes) |
| `--style` | `dot` | Mark style: `dot` (fixed-size marker, always visible) / `rect` (length-proportional) |
| `--dot-size` | `40` | Round marker area in points² (dot style only) |
| `--call-dist` | `30` | Telomere call: look within this many kb of each end |
| `--call-min` | `50` | Telomere call: forward+reverse in that end region ≥ this → capped |
| `--no-call` | off | Disable telomere calling (no end markers, no HTML report) |
| `--width` / `--height` | auto | Figure size (px) |
| `--dpi` | `200` | Raster (PNG) resolution |
| `--format` | `pdf` | `pdf` / `png` / `svg`, comma-separated |

## Development

```bash
pytest
```
