#!/usr/bin/env bash
#
# mosdepth_lowmapq_blacklist.sh — flag regions dominated by multi-mappers.
#
# Given paired mosdepth regions.bed.gz outputs, one unfiltered by MAPQ and
# one MAPQ-filtered, emit a BED of genomic windows where the fraction of
# reads passing the MAPQ filter is low enough that the locus is effectively
# "unmapped" for any caller that only trusts unique placements. These are
# prime candidate blacklist intervals for svaba.
#
# The input pair must come from the SAME BAM with the SAME --by, differing
# only in --mapq. Typical generation:
#
#   mosdepth --by 1000 --fast-mode --mapq  0 --threads 4 unf  SAMPLE.bam
#   mosdepth --by 1000 --fast-mode --mapq 30 --threads 4 filt SAMPLE.bam
#
# mosdepth writes half-open 0-based BED with 4 columns:
#   chrom  start  end  mean_depth
#
# This script joins the two files row-by-row (they are bin-aligned), computes
#   ratio = filt_mean_depth / unf_mean_depth
# and flags windows that pass BOTH of:
#   unf_mean_depth >= --min-cov         (skip near-empty bins where 0/0 dominates)
#   ratio          <  --max-ratio       (most reads are multi-mappers)
#
# Flagged per-bin rows are then bedtools-merged with --merge-gap into
# contiguous intervals, carrying the worst ratio and peak coverage forward
# so the merged output has useful label columns.
#
# The merged BED plugs straight into `combine_blacklists.sh` alongside your
# other component BEDs (tracks/hg38.blacklist.sorted.bed, high_runtime.bed,
# rmsk.simple_repeat.bed, etc.).
#
# Usage:
#   mosdepth_lowmapq_blacklist.sh [OPTIONS] UNF.regions.bed.gz FILT.regions.bed.gz
#
# Options:
#   --min-cov N         minimum unfiltered mean depth for a bin to be
#                       considered (default 10). Bins below this get skipped
#                       entirely — you do not want a 0/1 accidental flag on
#                       an empty region.
#   --max-ratio R       flag if filt/unf < R (default 0.5). Lower = stricter.
#                       Typical blacklist tuning range is 0.3 – 0.7.
#   --merge-gap N       bedtools merge -d value (default 1000). Flagged bins
#                       within N bp of each other get coalesced into one
#                       interval.
#   --out PREFIX        output prefix (default ./lowmapq).
#   --dump-per-region   also emit ${PREFIX}.per_region.tsv.gz with every
#                       considered bin (chr start end unf filt ratio). Off
#                       by default — on a WGS BAM this file is ~coverage-
#                       file-sized.
#   -h | --help         show this help
#
# Outputs (always):
#   ${PREFIX}.flagged.bed         one row per flagged bin, cols:
#                                   chr  start  end  ratio  unf_cov_rounded
#   ${PREFIX}.flagged.merged.bed  bedtools-merged intervals, cols:
#                                   chr  start  end  min_ratio  max_cov
#                                 (min_ratio = worst ratio inside the merged
#                                  block; max_cov = peak unfiltered coverage)
#   ${PREFIX}.summary.txt         coverage + flag tallies
#
# Outputs (optional):
#   ${PREFIX}.per_region.tsv.gz   with --dump-per-region
#
# Requires: gzip, awk, sort, bedtools on PATH.

set -euo pipefail

MIN_COV=10
MAX_RATIO=0.5
MERGE_GAP=1000
OUT_PREFIX=lowmapq
DUMP_PER=0

usage() {
  cat <<'EOF'
mosdepth_lowmapq_blacklist.sh — flag regions dominated by multi-mappers.

Usage:
  mosdepth_lowmapq_blacklist.sh [OPTIONS] UNF.regions.bed.gz FILT.regions.bed.gz

The two inputs must come from the SAME BAM with the SAME --by on mosdepth,
differing only in --mapq. Example generation:
  mosdepth --by 1000 --fast-mode --mapq  0 --threads 4 unf  sample.bam
  mosdepth --by 1000 --fast-mode --mapq 30 --threads 4 filt sample.bam

Options:
  --min-cov N         minimum unf mean-depth for a bin to be considered (default 10)
  --max-ratio R       flag if filt/unf < R (default 0.5); 0.3 stricter, 0.7 looser
  --merge-gap N       bedtools merge -d value (default 1000 bp)
  --out PREFIX        output prefix (default lowmapq)
  --dump-per-region   also write ${PREFIX}.per_region.tsv.gz
  -h, --help          show this message

Outputs:
  ${PREFIX}.flagged.bed         flagged bins:  chr start end ratio unf_cov
  ${PREFIX}.flagged.merged.bed  merged:        chr start end min_ratio max_cov
  ${PREFIX}.summary.txt         coverage + flag tallies

See the header comment in this script for the full rationale / column docs.
EOF
  exit 0
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --min-cov)         MIN_COV=$2; shift 2 ;;
    --max-ratio)       MAX_RATIO=$2; shift 2 ;;
    --merge-gap)       MERGE_GAP=$2; shift 2 ;;
    --out)             OUT_PREFIX=$2; shift 2 ;;
    --dump-per-region) DUMP_PER=1; shift ;;
    -h|--help)         usage ;;
    --)                shift; break ;;
    -*)                echo "unknown option: $1" >&2; exit 2 ;;
    *)                 break ;;
  esac
done

if [[ $# -ne 2 ]]; then
  echo "usage: $0 [opts] UNF.regions.bed.gz FILT.regions.bed.gz" >&2
  echo "       (try -h for full docs)" >&2
  exit 2
fi
UNF=$1
FILT=$2

for t in gzip awk sort bedtools; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }
done
for f in "$UNF" "$FILT"; do
  [[ -f $f ]] || { echo "ERROR: file not found: $f" >&2; exit 1; }
done

# Make sure the output directory exists if the user passed a pathy prefix.
OUT_DIR=$(dirname -- "$OUT_PREFIX")
[[ -d "$OUT_DIR" ]] || mkdir -p "$OUT_DIR"

echo "mosdepth_lowmapq_blacklist:"
echo "  unf=$UNF"
echo "  filt=$FILT"
echo "  min_cov=$MIN_COV  max_ratio=$MAX_RATIO  merge_gap=$MERGE_GAP"
echo "  out_prefix=$OUT_PREFIX  dump_per_region=$DUMP_PER"

# ---- 1. Join and filter -----------------------------------------------------
# Both mosdepth files were emitted from the same --by on the same reference,
# so their row-1..row-N bins are identical. We paste them and sanity-check
# the coords on every row — a mismatch usually means the user accidentally
# paired files from different BAMs or different --by values, which makes the
# ratio meaningless.
flagged="${OUT_PREFIX}.flagged.bed"
tmp_flag=$(mktemp)
tmp_per=$(mktemp)
trap 'rm -f "$tmp_flag" "$tmp_per" "${tmp_per}.gz"' EXIT

paste <(gzip -dc -- "$UNF") <(gzip -dc -- "$FILT") \
  | awk -F'\t' -v OFS='\t' \
        -v mc="$MIN_COV" -v mr="$MAX_RATIO" \
        -v flag="$tmp_flag" -v per="$tmp_per" -v dump="$DUMP_PER" '
    BEGIN {
      total       = 0
      considered  = 0
      flagged_cnt = 0
    }
    {
      # cols 1-4: unf (chrom start end meanDepth)
      # cols 5-8: filt (chrom start end meanDepth)
      if ($1 != $5 || $2 != $6 || $3 != $7) {
        printf "ERROR: row %d coordinate mismatch: %s:%s-%s vs %s:%s-%s\n" \
               "  likely --by mismatch or different BAMs.\n",
               NR, $1, $2, $3, $5, $6, $7 > "/dev/stderr"
        exit 1
      }
      total++
      u = $4 + 0
      f = $8 + 0
      if (u < mc) next
      considered++
      r = (u > 0 ? f / u : 0)
      if (dump == 1) {
        printf "%s\t%s\t%s\t%.3f\t%.3f\t%.4f\n", $1, $2, $3, u, f, r > per
      }
      if (r < mr) {
        # Col 4 = ratio (merge-friendly numeric), col 5 = rounded unf cov
        # so a `bedtools merge -c 4,5 -o min,max` gives "worst ratio in
        # block / peak coverage in block" — both useful for downstream
        # inspection / QC.
        printf "%s\t%s\t%s\t%.4f\t%d\n", $1, $2, $3, r, int(u + 0.5) > flag
        flagged_cnt++
      }
    }
    END {
      printf "scan: rows_total=%d  rows_considered=%d  rows_flagged=%d\n",
             total, considered, flagged_cnt > "/dev/stderr"
    }
  '

# ---- 2. Sort + write per-bin flagged BED ------------------------------------
# bedtools merge requires a sorted BED. mosdepth output is usually already in
# BAM-SQ order (which is often `chr1,chr2,...,chrX,chrY,chrM`), but we don't
# rely on that — re-sort with LC_ALL=C for deterministic output.
LC_ALL=C sort -k1,1 -k2,2n "$tmp_flag" > "$flagged"
n_flag=$(wc -l < "$flagged" | tr -d ' ')
echo "[1/3] wrote $flagged  ($n_flag flagged bins)"

# ---- 3. Merge ---------------------------------------------------------------
merged="${OUT_PREFIX}.flagged.merged.bed"
if [[ $n_flag -gt 0 ]]; then
  bedtools merge -i "$flagged" -d "$MERGE_GAP" \
    -c 4,5 -o min,max > "$merged"
else
  : > "$merged"
fi
n_merged=$(wc -l < "$merged" | tr -d ' ')
echo "[2/3] wrote $merged  ($n_merged merged intervals)"

# ---- 4. Summary -------------------------------------------------------------
summary="${OUT_PREFIX}.summary.txt"
flag_bp=$(awk   '{s+=$3-$2} END {print s+0}' "$flagged")
merged_bp=$(awk '{s+=$3-$2} END {print s+0}' "$merged")
{
  echo "mosdepth_lowmapq_blacklist  ($(date -u +%Y-%m-%dT%H:%M:%SZ))"
  echo "inputs:"
  echo "  unf:        $UNF"
  echo "  filt:       $FILT"
  echo "parameters:"
  echo "  min_cov:    $MIN_COV"
  echo "  max_ratio:  $MAX_RATIO"
  echo "  merge_gap:  $MERGE_GAP"
  echo ""
  printf "flagged bins:      %8d  (%12d bp,  %7.1f Mb)\n" \
         "$n_flag"   "$flag_bp"   "$(awk -v x="$flag_bp"   'BEGIN{printf "%.1f", x/1e6}')"
  printf "merged intervals:  %8d  (%12d bp,  %7.1f Mb)\n" \
         "$n_merged" "$merged_bp" "$(awk -v x="$merged_bp" 'BEGIN{printf "%.1f", x/1e6}')"
  echo ""
  echo "Top 20 widest merged intervals:"
  sort -k1,1 -k2,2n "$merged" \
    | awk 'BEGIN{OFS="\t"} {print $1, $2, $3, $3-$2, $4, $5}' \
    | sort -k4,4nr \
    | head -20 \
    | awk 'BEGIN{OFS="\t"; print "#chrom","start","end","width","min_ratio","max_cov"} {print}'
} > "$summary"
echo "[3/3] wrote $summary"

# ---- 5. Optional per-region dump -------------------------------------------
if [[ $DUMP_PER -eq 1 ]]; then
  out="${OUT_PREFIX}.per_region.tsv.gz"
  # Prepend header, then gzip. Keep it streaming so we don't materialize the
  # full file uncompressed.
  { printf "#chrom\tstart\tend\tunf_cov\tfilt_cov\tratio\n"; cat "$tmp_per"; } \
    | gzip -c > "$out"
  echo "     wrote $out"
fi

echo
cat "$summary"
