#!/usr/bin/env bash
#
# run_svaba_on_panel.sh -- run svaba (tumor/normal) on every simulated tumor BAM
# in a panel folder built by build_sim_panel.sh, then postprocess (+ optional
# benchmark vs the panel's ground truth). Mirrors the user's preferred invocation:
#
#   time svaba run -t <sim_tumor> -n <matched_normal> -G <ref> \
#        --blacklist hg38.combined_blacklist.bed -p 20 -a <id> \
#        --annotation hg38.rmsk.bed.gz --annotation hg38.segdups.bed.gz
#
# Each tumor BAM gets its own svaba_runs/<tag>/ output dir.
#
# Usage:
#   ./run_svaba_on_panel.sh [PANELDIR]        # default: newest /Volumes/wala24T/sim/wgs_*
#   NORMAL=... THREADS=20 ./run_svaba_on_panel.sh /Volumes/wala24T/sim/wgs_..._seed42
#   DRYRUN=1 ./run_svaba_on_panel.sh          # print the commands, run nothing
#
set -euo pipefail

# ---- config (override via env) ---------------------------------------------
: "${SVABA:=$HOME/git/svaba/build/svaba}"
: "${REF:=$HOME/ref_genome/Homo_sapiens_assembly38.fasta}"
: "${BLACKLIST:=$HOME/git/svaba/tracks/hg38.combined_blacklist.bed}"
: "${ANNO_RMSK:=$HOME/git/svaba/tracks/hg38.rmsk.bed.gz}"
: "${ANNO_SEGDUP:=$HOME/git/svaba/tracks/hg38.segdups.bed.gz}"
: "${THREADS:=20}"
: "${POSTPROCESS:=1}"           # run `svaba postprocess` after each (bps mode only)
: "${BENCHMARK:=1}"             # run sim/benchmark.py after each, write a summary
: "${BAM_GLOB:=tumor.*x*.bam}"  # which BAMs in the panel to run on (impure + pure)
: "${DRYRUN:=0}"

# Output mode = what we score against truth:
#   bps  (default) -- svaba2: postprocess -> <tag>.bps.sorted.dedup.txt.gz, scored
#                     with --pass-only --somatic-only.
#   vcf            -- older svaba (e.g. SVABA=~/git/svaba1/build/svaba): score the
#                     somatic SV + indel VCFs that `svaba run` writes directly
#                     (<tag>.svaba.somatic.{sv,indel}.vcf). No postprocess.
: "${MODE:=bps}"

# Pass --annotation to `svaba run`? auto (default) = on iff THIS binary supports
# it (older 2.0.x builds predate the flag and abort on it -- bps-mode comparisons
# of those builds need it off). Force with ANNOTATE=1 / ANNOTATE=0.
: "${ANNOTATE:=auto}"

# MATCHED NORMAL for the T/N pair. This is the user's SECOND WGS normal from the
# SAME individual as the tumor's contamination (build_sim_panel.sh's
# CONTAM_NORMAL_BAM = NP_WGS_Normal_DNA). Same person => shared germline (so
# germline SVs are correctly NOT somatic), but DIFFERENT reads from the
# contamination => no identical-read artifact. Only the simulated SVs are
# tumor-only/somatic. (Must be the same individual; a different person's normal
# would call that person's germline SVs as somatic.)
: "${NORMAL:=$HOME/Downloads/svaba_compare/blood.recal.bam}"

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
log(){ printf '[svaba-panel %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# ---- resolve panel dir -----------------------------------------------------
PANELDIR="${1:-}"
if [ -z "$PANELDIR" ]; then
  PANELDIR=$(ls -dt /Volumes/wala24T/sim/wgs_* 2>/dev/null | head -1 || true)
  [ -n "$PANELDIR" ] && log "no PANELDIR given -> newest: $PANELDIR"
fi
[ -n "$PANELDIR" ] && [ -d "$PANELDIR" ] || { echo "ERROR: panel dir not found (pass it as arg 1)"; exit 1; }
[ -f "$PANELDIR/truth.bedpe" ] || log "WARN: no truth.bedpe in $PANELDIR (benchmark will be skipped)"

# ---- preflight -------------------------------------------------------------
[ -x "$SVABA" ] || { echo "ERROR: svaba not executable: $SVABA"; exit 1; }
[ -f "$REF.fai" ] && [ -f "$REF.bwt" ] || { echo "ERROR: reference not indexed: $REF"; exit 1; }
[ -f "$NORMAL" ] || { echo "ERROR: matched normal not found: $NORMAL"; exit 1; }
[ -f "$NORMAL.bai" ] || log "WARN: matched normal has no .bai"
for b in "$BLACKLIST" "$ANNO_RMSK" "$ANNO_SEGDUP"; do
  [ -f "$b" ] || log "WARN: missing track $b (svaba will run without it if you drop the flag)"
done

shopt -s nullglob
BAMS=( "$PANELDIR"/$BAM_GLOB )
shopt -u nullglob
[ "${#BAMS[@]}" -gt 0 ] || { echo "ERROR: no BAMs matching '$BAM_GLOB' in $PANELDIR"; exit 1; }

# svaba build provenance: short git commit of the svaba source being run, so each
# commit's results land in their OWN folder (svaba_runs_<hash>) for tracking how
# accuracy/compute change across commits over time. Override with GITHASH=...
SVABA_REPO=$(cd "$(dirname "$SVABA")/.." 2>/dev/null && pwd || echo "")
if [ -z "${GITHASH:-}" ]; then
  if [ -n "$SVABA_REPO" ] && git -C "$SVABA_REPO" rev-parse --git-dir >/dev/null 2>&1; then
    GITHASH=$(git -C "$SVABA_REPO" rev-parse --short=8 HEAD)
    git -C "$SVABA_REPO" diff --quiet 2>/dev/null || GITHASH="${GITHASH}-dirty"
  else GITHASH="nogit"; fi
fi
SVABA_VER=$("$SVABA" --version 2>&1 | head -1 || true)
TIMEFORMAT='%R %U %S'              # bash `time` -> "real user sys" (seconds)
BATCH_START=$(date +%s)

# mode-dependent setup. --annotation is gated by ANNOTATE + binary support so a
# 2.0.x build (which aborts on the unknown flag) runs fine in bps mode.
case "$MODE" in
  bps)
    use_anno=0
    if [ "$ANNOTATE" = 1 ]; then use_anno=1
    elif [ "$ANNOTATE" = auto ]; then
      if "$SVABA" run --help 2>&1 | grep -q -- "--annotation"; then use_anno=1
      else log "NOTE: $SVABA has no --annotation (older build); running bps mode without it"; fi
    fi
    if [ "$use_anno" = 1 ]; then ANNO_ARGS=(--annotation "$ANNO_RMSK" --annotation "$ANNO_SEGDUP"); else ANNO_ARGS=(); fi ;;
  vcf) ANNO_ARGS=() ; POSTPROCESS=0 ;;
  *)   echo "ERROR: MODE must be 'bps' or 'vcf' (got '$MODE')"; exit 1 ;;
esac
anno_args(){ [ "${#ANNO_ARGS[@]}" -gt 0 ] && printf '%s ' "${ANNO_ARGS[@]}"; }
# combine svaba1's somatic SV + indel VCFs (plain or .gz) into one file for benchmark
combine_somatic_vcfs(){ local od="$1" tag="$2" out="$3"; : > "$out"; local v found=0
  shopt -s nullglob
  for v in "$od/$tag".svaba.somatic.sv.vcf "$od/$tag".svaba.somatic.sv.vcf.gz \
           "$od/$tag".svaba.somatic.indel.vcf "$od/$tag".svaba.somatic.indel.vcf.gz; do
    [ -f "$v" ] || continue; found=1
    case "$v" in *.gz) gzip -dc "$v";; *) cat "$v";; esac >> "$out"
  done
  shopt -u nullglob
  [ "$found" = 1 ]
}

RUNROOT="$PANELDIR/svaba_runs_${GITHASH}"
SUMMARY="$RUNROOT/benchmark_summary.tsv"
log "panel: $PANELDIR"
log "svaba: $SVABA   commit: $GITHASH   mode: $MODE   ($SVABA_VER)"
log "normal (-n): $NORMAL   threads: $THREADS   postprocess=$POSTPROCESS benchmark=$BENCHMARK"
log "tumor BAMs (${#BAMS[@]}): $(for b in "${BAMS[@]}"; do basename "$b"; done | tr '\n' ' ')"
log "output -> $RUNROOT/"

if [ "$DRYRUN" = "1" ]; then
  for bam in "${BAMS[@]}"; do tag=$(basename "${bam%.bam}")
    echo "time $SVABA run -t $bam -n $NORMAL -G $REF --blacklist $BLACKLIST -p $THREADS -a $tag $(anno_args)"
  done
  log "DRYRUN=1 (mode=$MODE) -> nothing executed."; exit 0
fi

mkdir -p "$RUNROOT"
printf "tumor_bam\ttruth\tdetected\trecall\tcalls\tmatching\tprecision\twall_s\tcpu_s\n" > "$SUMMARY"

# ---- run svaba on each tumor BAM -------------------------------------------
for bam in "${BAMS[@]}"; do
  tag=$(basename "${bam%.bam}")
  od="$RUNROOT/$tag"; mkdir -p "$od"
  primary=$([ "$MODE" = vcf ] && echo "$od/$tag.svaba.somatic.sv.vcf" || echo "$od/$tag.bps.txt.gz")
  if [ -f "$primary" ]; then
    log "$tag: svaba output exists -> skip run"
  else
    log "svaba run on $tag ($(basename "$bam")) [mode=$MODE] ..."
    # capture real/user/sys for this run into <tag>.time.txt (svaba's own output
    # goes to the log); the `time` keyword's stderr is what lands in time.txt.
    { time ( cd "$od" && "$SVABA" run -t "$bam" -n "$NORMAL" -G "$REF" \
        --blacklist "$BLACKLIST" -p "$THREADS" -a "$tag" \
        ${ANNO_ARGS[@]+"${ANNO_ARGS[@]}"} \
        > "$od/$tag.svaba.log" 2>&1 ) ; } 2> "$od/$tag.time.txt" \
      || log "  WARN: svaba run nonzero exit for $tag (see $tag.svaba.log)"
    log "  done ($(cat "$od/$tag.time.txt" 2>/dev/null | awk '{printf "wall %ss cpu %.0fs", $1, $2+$3}'))"
  fi

  if [ "$POSTPROCESS" = "1" ]; then
    if [ -f "$od/$tag.bps.sorted.dedup.txt.gz" ]; then
      log "  postprocess output exists -> skip"
    else
      log "  postprocess $tag ..."
      ( cd "$od" && "$SVABA" postprocess -i "$tag" -t "$THREADS" ) >> "$od/$tag.svaba.log" 2>&1 || \
        log "  WARN: postprocess failed (see log)"
    fi
  fi

  # per-run compute (real wall, user+sys cpu) from <tag>.time.txt
  wall="NA"; cpu="NA"
  if [ -s "$od/$tag.time.txt" ]; then
    read -r rrun urun srun < "$od/$tag.time.txt" 2>/dev/null || true
    wall="${rrun:-NA}"
    cpu=$(awk -v u="${urun:-0}" -v s="${srun:-0}" 'BEGIN{printf "%.1f",u+s}')
  fi

  # accuracy vs truth (optional). bps mode -> dedup bps (--somatic-only); vcf mode
  # -> combined somatic SV+indel VCFs (already the somatic set, so no --somatic-only).
  t="NA"; d="NA"; r="NA"; c="NA"; m="NA"; pr="NA"; calls=""; bflags=""
  if [ "$MODE" = "vcf" ]; then
    if combine_somatic_vcfs "$od" "$tag" "$od/$tag.somatic.combined.vcf"; then
      calls="$od/$tag.somatic.combined.vcf"; bflags="--format vcf --pass-only"
    else log "  WARN: no somatic VCFs found for $tag (expected $tag.svaba.somatic.{sv,indel}.vcf)"; fi
  else
    [ -f "$od/$tag.bps.sorted.dedup.txt.gz" ] && { calls="$od/$tag.bps.sorted.dedup.txt.gz"; bflags="--pass-only --somatic-only"; }
  fi
  if [ "$BENCHMARK" = "1" ] && [ -f "$PANELDIR/truth.bedpe" ] && [ -n "$calls" ]; then
    log "  benchmarking $tag vs truth ($MODE) ..."
    out=$(python3 "$HERE/benchmark.py" --truth "$PANELDIR/truth.bedpe" \
           --calls "$calls" $bflags --events-out "$od/$tag.events.tsv" 2>/dev/null || true)
    echo "$out" > "$od/$tag.benchmark.txt"
    ov=$(echo "$out" | grep -m1 "truth=[0-9]"); cv=$(echo "$out" | grep -m1 "calls=[0-9]")
    t=$(echo "$ov"  | sed -n 's/.*truth=\([0-9]*\).*/\1/p');  t=${t:-NA}
    d=$(echo "$ov"  | sed -n 's/.*detected=\([0-9]*\).*/\1/p'); d=${d:-NA}
    r=$(echo "$ov"  | sed -n 's/.*recall=\([0-9.]*%\).*/\1/p'); r=${r:-NA}
    c=$(echo "$cv"  | sed -n 's/.*calls=\([0-9]*\).*/\1/p');  c=${c:-NA}
    m=$(echo "$cv"  | sed -n 's/.*matching-truth=\([0-9]*\).*/\1/p'); m=${m:-NA}
    pr=$(echo "$cv" | sed -n 's/.*precision=\([0-9.]*%\).*/\1/p'); pr=${pr:-NA}
  fi
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$tag" "$t" "$d" "$r" "$c" "$m" "$pr" "$wall" "$cpu" >> "$SUMMARY"
done

# ---- totals + run metadata (for the cross-commit comparison tool) ----------
shopt -s nullglob; TF=( "$RUNROOT"/*/*.time.txt ); shopt -u nullglob
TOT_WALL=0; TOT_CPU=0
if [ "${#TF[@]}" -gt 0 ]; then
  # note the trailing \n + `|| true`: a newline-less `read` returns non-zero,
  # which under `set -e` would otherwise abort before run_meta.json is written.
  read -r TOT_WALL TOT_CPU < <(awk '{w+=$1;c+=$2+$3} END{printf "%.1f %.1f\n",w,c}' "${TF[@]}") || true
fi
BATCH_WALL=$(( $(date +%s) - BATCH_START ))
cat > "$RUNROOT/run_meta.json" <<JSON
{
  "svaba_git": "$GITHASH",
  "svaba_version": "$(printf '%s' "$SVABA_VER" | sed 's/\\/\\\\/g; s/"/\\"/g')",
  "date_utc": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "host": "$(hostname)",
  "panel": "$PANELDIR",
  "normal": "$NORMAL",
  "threads": $THREADS,
  "n_bams": ${#BAMS[@]},
  "total_svaba_wall_s": $TOT_WALL,
  "total_svaba_cpu_s": $TOT_CPU,
  "batch_wall_s": $BATCH_WALL
}
JSON

log "ALL DONE. svaba outputs in $RUNROOT/  (run_meta.json + benchmark_summary.tsv written)"
log "TOTAL svaba compute over ${#BAMS[@]} runs: wall=${TOT_WALL}s  cpu=${TOT_CPU}s  (batch wall incl. postprocess/benchmark=${BATCH_WALL}s)"
if [ -f "$SUMMARY" ]; then echo; column -t -s$'\t' "$SUMMARY" >&2 || cat "$SUMMARY" >&2; fi
log "Interactive scoring: open docs/sim_benchmark.html, load $PANELDIR/truth.bedpe + a <tag>.bps.sorted.dedup.txt.gz"
log "(Each svaba commit writes its own svaba_runs_<hash>/ + run_meta.json -> ready for the cross-commit comparison tool.)"
