#!/usr/bin/env bash
# thread_scaling.sh — measure svaba run wall / CPU time vs thread count.
#
# Runs `svaba run` (tumor/normal, one region) once per thread count and
# appends a row per run to $OUT/thread_scaling.tsv:
#   threads wall_s user_s sys_s cpu_s max_rss_kb pct_cpu n_bps exit
# Plot with: Rscript opt/thread_scaling_plot.R $OUT/thread_scaling.tsv opt/thread_scaling/
#
# Env overrides: REGION, THREAD_LIST, OUT, TUMOR, NORMAL, REF, BLACKLIST, SVABA
set -uo pipefail

SVABA_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
SVABA=${SVABA:-$SVABA_DIR/build-jemalloc/svaba}
REF=${REF:-$HOME/ref/Homo_sapiens_assembly38.fasta}
BLACKLIST=${BLACKLIST:-$SVABA_DIR/tracks/hg38.combined_blacklist.bed}
TUMOR=${TUMOR:-/mnt/wgs/SG.wgs.UCLA.2025.01.tumor_cleaned.recal.bam}
NORMAL=${NORMAL:-/mnt/wgs/blood.ss.cleaned.bam}
REGION=${REGION:-chr17}
THREAD_LIST=${THREAD_LIST:-"4 6 8 10 12 14 16 18 20 22 24 26 28 30 32"}
OUT=${OUT:-/mnt/wgs/svaba_bench}

mkdir -p "$OUT"; cd "$OUT"
TSV=$OUT/thread_scaling.tsv
[[ -f $TSV ]] || printf "threads\twall_s\tuser_s\tsys_s\tcpu_s\tmax_rss_kb\tpct_cpu\tn_bps\texit\n" > "$TSV"

# Record the environment once so the plot/notes can cite it.
{
  echo "date: $(date -u +%FT%TZ)"
  echo "host: $(hostname)"
  echo "svaba: $SVABA ($(cd "$SVABA_DIR" && git rev-parse --short HEAD))"
  echo "region: $REGION"
  echo "tumor: $TUMOR"; echo "normal: $NORMAL"
  echo "allocator: jemalloc via svaba_jemalloc"
  lscpu | grep -E 'Model name|^CPU\(s\)|Thread|Core|Socket|L3|NUMA node\(s\)'
  free -g | head -2
} > "$OUT/environment.txt"

# Warm the page cache: reference index + the region's reads, so run #1
# isn't penalized by cold disk reads relative to later runs.
echo "[$(date)] warming page cache"
for f in "$REF" "$REF.bwt" "$REF.sa" "$REF.pac"; do cat "$f" > /dev/null; done
samtools view -c -@4 "$TUMOR" "$REGION" > /dev/null
samtools view -c -@4 "$NORMAL" "$REGION" > /dev/null

for p in $THREAD_LIST; do
  id="scale_p$(printf %02d "$p")"
  if awk -v p="$p" 'NR>1 && $1==p && $9==0 {f=1} END {exit !f}' "$TSV"; then
    echo "[$(date)] p=$p already done, skipping"; continue
  fi
  echo "[$(date)] starting p=$p"
  rm -f "$id".*
  /usr/bin/time -o "$id.time" -f "%e\t%U\t%S\t%M\t%P" \
    env SVABA="$SVABA" \
        MALLOC_CONF="background_thread:true,narenas:$p,dirty_decay_ms:10000" \
    "$SVABA_DIR/svaba_jemalloc" run \
      -t "$TUMOR" -n "$NORMAL" -G "$REF" --blacklist "$BLACKLIST" \
      -k "$REGION" -p "$p" -a "$id" > "$id.stdout" 2>&1
  rc=$?
  IFS=$'\t' read -r wall user sys rss pct < <(tail -1 "$id.time")
  nbps=$( (zcat "$id.bps.txt.gz" 2>/dev/null || true) | grep -vc '^#' )
  cpu=$(awk -v u="$user" -v s="$sys" 'BEGIN{printf "%.2f", u+s}')
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$p" "$wall" "$user" "$sys" "$cpu" "$rss" "${pct%\%}" "$nbps" "$rc" >> "$TSV"
  echo "[$(date)] p=$p wall=${wall}s cpu=${cpu}s rss=$((rss/1024))MB exit=$rc"
done
echo "[$(date)] all done -> $TSV"
