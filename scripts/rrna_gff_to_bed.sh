#!/usr/bin/env bash
# rrna_gff_to_bed.sh — turn a barrnap / pybarrnap rRNA GFF into a teloviz feature
# BED (rDNA/NOR arrays). Run this ONCE before the main teloviz run:
#
#   pybarrnap  ->  rRNA.gff
#      │ (this script)
#      ▼
#   features.bed  ->  teloviz ... --rDNA features.bed
#
# This is a preprocessing companion, NOT part of the `teloviz` command itself:
# teloviz only ever reads BED (one parser, one contract). The biological curation
# lives HERE, out in the open.
#
# 45S (the NOR):  merge 18S/5.8S/28S hits, then keep an interval ONLY if all
#   three components are present. A complete 45S unit is 18S-ITS1-5.8S-ITS2-28S,
#   so "are the three there?" is a BINARY test — no count threshold to tune (and
#   a count cutoff would wrongly drop a cleanly-detected single unit, which is 3
#   hits). This is what tells a real NOR from a stray 18S/28S fragment.
#
# 5S (optional, --with-5s):  a single gene with no internal composition, so the
#   only handle is copy count. Off by default because the 5S array often sits
#   mid-chromosome, irrelevant to end-cap QC. Enable to add it as a lane.
#
# Requires: bedtools, awk, grep, sort (all standard).
#
# The output BED's chromosome names come from the GFF and MUST match the names in
# your tidk windows TSV and .fai, or teloviz drops every feature with a
# "dropped N features" warning. Rename one side first if they differ
# (e.g. chr1 vs CM061962.1).

set -euo pipefail

# ---- defaults ---------------------------------------------------------------
MERGE=50000        # bp: bridge tandem rRNA copies into one array. 45S units sit
                   # ~31 kb apart (IGS); 35 kb-500 kb all give the same intervals
                   # (a plateau), so 50000 is a representative value, not magic.
WITH_5S=0          # include the 5S array lane? (off by default)
MIN_5S=10          # with --with-5s: keep 5S clusters with >= this many copies
                   # (a real array has many copies; scattered singletons are
                   # 1-2, so the count distribution is usually clearly bimodal)
DRY_RUN=0

usage() {
  cat <<'EOF'
Usage: rrna_gff_to_bed.sh [options] <input.gff> <output.bed>
       rrna_gff_to_bed.sh --dry-run [options] <input.gff>

Convert a barrnap/pybarrnap rRNA GFF into a teloviz feature BED.
  rdna_45S : 18S/5.8S/28S merged into arrays, kept only where ALL THREE
             components are present (composition filter, no count threshold).
  rdna_5S  : (only with --with-5s) 5S merged into arrays, clusters below the
             copy-count cutoff dropped.

Options:
  --merge N     Merge distance for tandem copies, bp        [default 50000]
  --with-5s     Also emit the 5S array lane (off by default)
  --5s-min N    With --with-5s: keep 5S arrays with >= N copies [default 10]
  --dry-run     Show the 45S interval composition and the 5S count distribution
                (to sanity-check the cutoffs), then exit without writing a BED.
  -h, --help    This help.

Notes:
  * 45S needs NO threshold: a complete unit is exactly 18S+5.8S+28S, so the
    filter is binary. Do not "improve" it with a count cutoff — it would drop
    clean single units (chr17/chr32 are single units of 4 hits).
  * Sanity check: the rdna_45S features should land where you expect the NOR
    (often near the ends of acrocentric chromosomes).
EOF
}

# ---- arg parsing ------------------------------------------------------------
POS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --merge)   MERGE="$2"; shift 2 ;;
    --with-5s) WITH_5S=1;  shift ;;
    --5s-min)  MIN_5S="$2"; shift 2 ;;
    --dry-run) DRY_RUN=1;  shift ;;
    -h|--help) usage; exit 0 ;;
    --) shift; while [[ $# -gt 0 ]]; do POS+=("$1"); shift; done ;;
    -*) echo "rrna_gff_to_bed: unknown option: $1" >&2; usage >&2; exit 2 ;;
    *)  POS+=("$1"); shift ;;
  esac
done

if ! command -v bedtools >/dev/null 2>&1; then
  echo "rrna_gff_to_bed: bedtools not found on PATH (conda install -c bioconda bedtools)." >&2
  exit 3
fi

GFF="${POS[0]:-}"
if [[ -z "$GFF" || ! -f "$GFF" ]]; then
  echo "rrna_gff_to_bed: input GFF not found: '${GFF:-<none>}'" >&2
  usage >&2; exit 2
fi

# 45S merged intervals with copy count + the distinct component set. Reused by
# both the dry-run summary and the real output.
merged_45s() {
  grep -E "Name=(18S|28S|5_8S)_rRNA" "$GFF" \
    | awk 'BEGIN{OFS="\t"} !/^#/ {match($9,/Name=[^;]+/);
           print $1,$4-1,$5,substr($9,RSTART+5,RLENGTH-5)}' \
    | sort -k1,1 -k2,2n \
    | bedtools merge -d "$MERGE" -c 4,4 -o count,distinct
}

# ---- dry-run: show what the cutoffs are doing, write nothing -----------------
if [[ "$DRY_RUN" -eq 1 ]]; then
  echo "== 45S merged intervals (merge -d ${MERGE}): PASS = 18S+5.8S+28S all present =="
  merged_45s | awk 'BEGIN{OFS="\t"}
    {pass=($5~/18S_rRNA/ && $5~/28S_rRNA/ && $5~/5_8S_rRNA/) ? "PASS" : "drop";
     print pass,$1,$2,$3,"n="$4,$5}'
  echo
  echo "== 5S cluster sizes (merge -d ${MERGE}), copies per cluster, high -> low =="
  echo "   (only used with --with-5s; put --5s-min in the gap):"
  grep "Name=5S_rRNA" "$GFF" \
    | awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$4-1,$5}' | sort -k1,1 -k2,2n \
    | bedtools merge -d "$MERGE" -c 1 -o count \
    | awk '{print $4}' | sort -rn | tr '\n' ' '
  echo
  exit 0
fi

OUT="${POS[1]:-}"
if [[ -z "$OUT" ]]; then
  echo "rrna_gff_to_bed: output BED path required (or use --dry-run)." >&2
  usage >&2; exit 2
fi

tmp="$(mktemp)"
trap 'rm -f "$tmp"' EXIT

# 45S: keep only intervals containing all three components. (Pattern and action
# MUST stay on one line — a newline between them makes awk two separate rules.)
merged_45s \
  | awk 'BEGIN{OFS="\t"} $5~/18S_rRNA/ && $5~/28S_rRNA/ && $5~/5_8S_rRNA/ {print $1,$2,$3,"45S_n"$4,"rdna_45S"}' > "$tmp"

# 5S (optional): arrays only, scattered clusters below the cutoff dropped.
if [[ "$WITH_5S" -eq 1 ]]; then
  grep "Name=5S_rRNA" "$GFF" \
    | awk 'BEGIN{OFS="\t"} !/^#/ {print $1,$4-1,$5}' | sort -k1,1 -k2,2n \
    | bedtools merge -d "$MERGE" -c 1 -o count \
    | awk -v m="$MIN_5S" 'BEGIN{OFS="\t"} $4>=m {print $1,$2,$3,"5S_n"$4,"rdna_5S"}' >> "$tmp"
fi

sort -k1,1 -k2,2n "$tmp" -o "$OUT"

n45=$(awk '$5=="rdna_45S"' "$OUT" | wc -l | tr -d ' ')
n5=$(awk '$5=="rdna_5S"' "$OUT" | wc -l | tr -d ' ')
echo "rrna_gff_to_bed: wrote $OUT"
echo "  rdna_45S arrays (18S+5.8S+28S complete): $n45   (merge=${MERGE})"
if [[ "$WITH_5S" -eq 1 ]]; then
  echo "  rdna_5S  arrays: $n5    (merge=${MERGE}, 5s-min=${MIN_5S})"
else
  echo "  rdna_5S: skipped (pass --with-5s to include)"
fi
echo "  Next: teloviz <windows.tsv> --fai <genome.fai> --rDNA $OUT -o <prefix>"
echo "  (chromosome names in $OUT must match the tidk TSV / .fai.)"
