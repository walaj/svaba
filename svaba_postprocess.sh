#!/usr/bin/env bash
#
# svaba_postprocess.sh — one-stop post-processing of a svaba run.
#
# Replaces (and subsumes) the older sort_output.sh and
# sort_and_deduplicate_bps.sh. Given a svaba analysis ID, this script:
#
#   1. Merges per-thread output BAMs into one BAM per suffix.
#        ${ID}.thread*.${suffix}.bam -> ${ID}.${suffix}.bam
#
#   2. Runs `svaba postprocess` (the C++ subcommand) to coord-sort, stream-
#      dedup on exact (qname,flag) duplicates, stamp a @PG "svaba_postprocess"
#      line onto every final BAM header, and build a .bai index. Suffixes:
#        weird, corrected, discordant, contigs.
#      Final state per suffix is ${ID}.${suffix}.bam(.bai); no
#      .sorted / .deduped intermediate files remain on the user path.
#
#   3. Sorts, deduplicates, and PASS-filters ${ID}.bps.txt.gz. Emits:
#        ${ID}.bps.sorted.txt.gz            — full sorted list
#        ${ID}.bps.sorted.dedup.txt.gz      — best SV (max maxlod) per
#                                              (chr1,pos1,strand1,chr2,pos2,strand2)
#        ${ID}.bps.sorted.dedup.pass.txt.gz — dedup restricted to PASS rows
#
#   4. Filters ${ID}.r2c.txt.gz down to the contigs of PASS breakpoints,
#      writing ${ID}.r2c.pass.txt.gz for interactive use in bps_explorer.html
#      (the full r2c.txt.gz is typically too heavy to load).
#
#   5. Optionally (--split-by-source) demultiplexes the deduped BAMs by the
#      first 4 chars of each read's QNAME into per-source BAMs.
#
# Every step is idempotent — missing inputs cause that step to log a short
# note and continue, so rerunning this after regenerating one of the inputs
# Just Works. The --skip-* flags let you skip a stage explicitly when you
# only want to refresh one phase.
#
# Usage:
#   svaba_postprocess.sh [options] <ID>
#
# Options (flags take precedence over equivalent env vars):
#   -t, --threads N       total threads budget                 (default: 4)
#   -m, --mem STR         per-samtools-sort-thread memory      (default: 2G)
#       --sort-buffer STR per-chunk memory for the bps sorter  (default: 2G)
#       --split-by-source also split BAMs by QNAME source prefix
#       --input-dir DIR   look for inputs in DIR               (default: .)
#       --output-dir DIR  write outputs to DIR               (default: INPUT_DIR)
#       --svaba PATH      path to svaba binary (default: ~/git/svaba/build/svaba)
#       --keep-tmp        keep scratch dirs for debugging
#       --skip-bam        skip BAM merge/sort/dedup/index (steps 1 & 2)
#       --skip-bps        skip bps sort/dedup/PASS (step 3)
#       --skip-r2c        skip r2c PASS filter (step 4)
#       --skip-split      force-disable split-by-source
#   -h, --help            this message
#
# Environment variable fallbacks (used only if the matching flag is absent):
#   THREADS, MEM, BUFFER_SIZE (→ --sort-buffer), SVABA, SAM, INPUT_DIR,
#   OUTPUT_DIR, SPLIT_BY_SOURCE, KEEP_TMP
#
# Examples:
#   svaba_postprocess.sh tumor_vs_normal
#   svaba_postprocess.sh -t 16 -m 4G tumor_vs_normal
#   svaba_postprocess.sh --split-by-source --output-dir /tmp/run42 run42
#   svaba_postprocess.sh --skip-bam --skip-bps tumor_vs_normal   # only r2c refresh
#
# Notes:
#   * Step 3 needs GNU sort (or `gsort` on macOS via `brew install coreutils`).
#   * Steps 1+2 need samtools on PATH; step 2 needs a working svaba binary.

set -euo pipefail

# ----------------------------------------------------------- defaults ---
THREADS=${THREADS:-4}
MEM=${MEM:-2G}
SORT_BUFFER=${BUFFER_SIZE:-${SORT_BUFFER:-2G}}
SVABA=${SVABA:-~/git/svaba/build/svaba}
SAM=${SAM:-samtools}
SPLIT_BY_SOURCE=${SPLIT_BY_SOURCE:-0}
INPUT_DIR=${INPUT_DIR:-.}
OUTPUT_DIR=${OUTPUT_DIR:-}
KEEP_TMP=${KEEP_TMP:-0}
SKIP_BAM=0
SKIP_BPS=0
SKIP_R2C=0

print_usage() {
  sed -n '2,65p' "$0" | sed 's/^# \{0,1\}//'
}

# ---------------------------------------------------------- arg parse ---
while [[ $# -gt 0 ]]; do
  case "$1" in
    -t|--threads)      THREADS="$2"; shift 2 ;;
    -m|--mem)          MEM="$2"; shift 2 ;;
    --sort-buffer)     SORT_BUFFER="$2"; shift 2 ;;
    --split-by-source) SPLIT_BY_SOURCE=1; shift ;;
    --skip-split)      SPLIT_BY_SOURCE=0; shift ;;
    --input-dir)       INPUT_DIR="$2"; shift 2 ;;
    --output-dir)      OUTPUT_DIR="$2"; shift 2 ;;
    --svaba)           SVABA="$2"; shift 2 ;;
    --keep-tmp)        KEEP_TMP=1; shift ;;
    --skip-bam)        SKIP_BAM=1; shift ;;
    --skip-bps)        SKIP_BPS=1; shift ;;
    --skip-r2c)        SKIP_R2C=1; shift ;;
    -h|--help)         print_usage; exit 0 ;;
    --)                shift; break ;;
    -*)
      echo "svaba_postprocess.sh: unknown option '$1'" >&2
      echo "run 'svaba_postprocess.sh --help' for usage" >&2
      exit 2
      ;;
    *)                 break ;;
  esac
done

if [[ $# -lt 1 ]]; then
  print_usage
  exit 2
fi
ID=$1
shift

OUTPUT_DIR=${OUTPUT_DIR:-$INPUT_DIR}
if [[ ! -d "$INPUT_DIR" ]]; then
  echo "svaba_postprocess.sh: INPUT_DIR does not exist: $INPUT_DIR" >&2
  exit 3
fi
mkdir -p "$OUTPUT_DIR"

# Human-readable echo of the plan so logs tell you what was about to run.
echo "svaba_postprocess.sh: id=$ID threads=$THREADS mem=$MEM sort_buffer=$SORT_BUFFER"
echo "                    : input_dir=$INPUT_DIR output_dir=$OUTPUT_DIR"
echo "                    : split_by_source=$SPLIT_BY_SOURCE skip_bam=$SKIP_BAM skip_bps=$SKIP_BPS skip_r2c=$SKIP_R2C"

# Useful for absolute-path filenames in the steps below. Use trailing-slash-
# free form so we always concatenate `${DIR}/${name}`.
INPUT_DIR=${INPUT_DIR%/}
OUTPUT_DIR=${OUTPUT_DIR%/}

# ========================================================================
# STEP 1: Merge per-thread BAMs.
#
# svaba run with multiple threads may produce ${ID}.thread0.${suffix}.bam,
# ${ID}.thread1.${suffix}.bam, ... Coalesce them into a single
# ${ID}.${suffix}.bam (to match what step 2 expects). Single-file inputs
# are renamed; no-file inputs are silently skipped.
# ========================================================================
if [[ $SKIP_BAM -eq 0 ]]; then
  echo "[1/5] merging per-thread BAMs"
  for suffix in discordant weird corrected; do
    target="${INPUT_DIR}/${ID}.${suffix}.bam"
    shopt -s nullglob
    bam_files=("${INPUT_DIR}/${ID}".thread*."${suffix}".bam)
    shopt -u nullglob
    if [[ ${#bam_files[@]} -gt 1 ]]; then
      echo "   merge ${#bam_files[@]} thread BAMs for '$suffix' -> $target"
      "$SAM" merge -f -@ "$THREADS" "$target" "${bam_files[@]}"
      rm -f "${bam_files[@]}"
    elif [[ ${#bam_files[@]} -eq 1 ]]; then
      echo "   mv  ${bam_files[0]}  ->  $target"
      mv "${bam_files[0]}" "$target"
    fi
  done
else
  echo "[1/5] skipping BAM merge (--skip-bam)"
fi

# ========================================================================
# STEP 2: svaba postprocess — sort, stream-dedup, stamp @PG, build .bai.
#
# Heavy lifting lives in the C++ subcommand. We just pass it the thread /
# memory budget; svaba postprocess splits THREADS across suffixes and passes
# -m to samtools sort per-thread.
# ========================================================================
if [[ $SKIP_BAM -eq 0 ]]; then
  echo "[2/5] svaba postprocess: sort + dedup + @PG + BAI index (threads=$THREADS, mem=$MEM)"
  (cd "$INPUT_DIR" && "$SVABA" postprocess -i "$ID" -t "$THREADS" -m "$MEM")
else
  echo "[2/5] skipping svaba postprocess (--skip-bam)"
fi

# ========================================================================
# STEP 3: Sort, dedup, and PASS-filter the bps.txt.
#
# Three outputs: sorted, sorted+dedup (one row per unique breakpoint pair,
# keeping the row with highest maxlod), sorted+dedup+PASS. Preserves the
# leading "#chr1..." header row in every output. Uses GNU sort (or gsort
# on macOS) for external-memory sort; --parallel=THREADS, -S SORT_BUFFER.
#
# Column positions in bps.txt (hard-coded from BreakPoint::toFileString):
#   col 32 = confidence (PASS / LOWMAPQ / BLACKLIST / ...)
#   col 38 = maxlod     (numeric; "NA" sorts last with -g)
#
# Sort keys in order: chr1(V), pos1(n), strand1, chr2(V), pos2(n), strand2,
#                     maxlod(gr, descending).
# ========================================================================
BPS_GZ="${INPUT_DIR}/${ID}.bps.txt.gz"
BPS_TXT="${INPUT_DIR}/${ID}.bps.txt"
OUT_SORTED="${OUTPUT_DIR}/${ID}.bps.sorted.txt.gz"
OUT_DEDUP="${OUTPUT_DIR}/${ID}.bps.sorted.dedup.txt.gz"
OUT_PASS="${OUTPUT_DIR}/${ID}.bps.sorted.dedup.pass.txt.gz"
CONF_COL=32
MAXLOD_COL=38

pick_sort() {
  if command -v gsort >/dev/null 2>&1; then
    echo gsort; return 0
  fi
  if sort --version 2>/dev/null | head -n 1 | grep -qi 'gnu'; then
    echo sort; return 0
  fi
  return 1
}

if [[ $SKIP_BPS -eq 0 ]]; then
  if [[ -r "$BPS_GZ"  ]]; then IN_BPS="$BPS_GZ";  READER="gzip -dc"
  elif [[ -r "$BPS_TXT" ]]; then IN_BPS="$BPS_TXT"; READER="cat"
  else IN_BPS=""; fi

  if [[ -z "$IN_BPS" ]]; then
    echo "[3/5] skipping bps sort/dedup/PASS: $BPS_GZ (or .bps.txt) not found"
  elif ! SORT_BIN=$(pick_sort); then
    echo "[3/5] skipping bps sort/dedup/PASS: GNU sort not on PATH" >&2
    echo "      (macOS: 'brew install coreutils' for gsort)" >&2
  else
    echo "[3/5] sort + dedup + PASS-filter  $IN_BPS"
    echo "      sort=$SORT_BIN  buffer=$SORT_BUFFER  parallel=$THREADS"

    TMP_BPS=$(mktemp -d "${TMPDIR:-/tmp}/svaba_pp_bps.XXXXXXXX")
    # Single cleanup trap for the whole script; registered the first time
    # we need a tempdir. Registered here because step 3 creates the first
    # tempdir used.
    trap 'if [[ "${KEEP_TMP:-0}" == "1" ]]; then
            echo "svaba_postprocess.sh: leaving scratch at: $TMP_BPS${TMP_SPLIT:+ $TMP_SPLIT}" >&2
          else
            rm -rf "$TMP_BPS" "${TMP_SPLIT:-}"
          fi' EXIT INT TERM

    # Peek at the first line to decide whether to split off a header.
    FIRST_LINE=$($READER "$IN_BPS" | head -n 1 || true)
    HAS_HEADER=0
    [[ "${FIRST_LINE:0:1}" == "#" ]] && HAS_HEADER=1

    # --- 3a) sort ---
    {
      if [[ $HAS_HEADER -eq 1 ]]; then
        printf '%s\n' "$FIRST_LINE"
        $READER "$IN_BPS" | tail -n +2
      else
        $READER "$IN_BPS"
      fi | {
        if [[ $HAS_HEADER -eq 1 ]]; then
          IFS= read -r hdr
          printf '%s\n' "$hdr"
        fi
        LC_ALL=C "$SORT_BIN" \
          -t $'\t' \
          -k1,1V \
          -k2,2n \
          -k3,3 \
          -k4,4V \
          -k5,5n \
          -k6,6 \
          -k${MAXLOD_COL},${MAXLOD_COL}gr \
          -S "$SORT_BUFFER" \
          --parallel="$THREADS" \
          -T "$TMP_BPS"
      }
    } | gzip -c > "$OUT_SORTED"
    echo "      wrote $OUT_SORTED"

    # --- 3b) dedup by breakpoint pair (first row per key wins = highest maxlod) ---
    gzip -dc "$OUT_SORTED" | awk -F'\t' -v OFS='\t' '
      NR == 1 && $0 ~ /^#/ { print; next }
      { key = $1 FS $2 FS $3 FS $4 FS $5 FS $6
        if (!(key in seen)) { seen[key] = 1; print } }
    ' | gzip -c > "$OUT_DEDUP"
    echo "      wrote $OUT_DEDUP"

    # --- 3c) PASS-filter the dedup ---
    gzip -dc "$OUT_DEDUP" | awk -F'\t' -v OFS='\t' -v CCOL="$CONF_COL" '
      NR == 1 && $0 ~ /^#/ { print; next }
      $CCOL == "PASS" { print }
    ' | gzip -c > "$OUT_PASS"
    echo "      wrote $OUT_PASS"
  fi
else
  echo "[3/5] skipping bps sort/dedup/PASS (--skip-bps)"
fi

# ========================================================================
# STEP 4: Filter r2c.txt.gz to PASS contigs only.
#
# The structured r2c TSV emitted by `svaba run` contains one row per
# assembled contig + one row per r2c-aligned read, for every variant-
# bearing contig. At full genome scale this can be tens of millions of
# rows — too heavy for the interactive viewer. Here we resolve a PASS
# contig set from bps.txt.gz (col 32 == "PASS", col 30 == cname) and
# keep only r2c rows whose contig_name (col 2) is in that set. Both the
# contig row and its read rows survive together because they share col 2.
# ========================================================================
R2C_GZ="${INPUT_DIR}/${ID}.r2c.txt.gz"
R2C_OUT="${OUTPUT_DIR}/${ID}.r2c.pass.txt.gz"

if [[ $SKIP_R2C -eq 0 ]]; then
  if [[ ! -f "$BPS_GZ" ]]; then
    echo "[4/5] skipping r2c PASS filter: $BPS_GZ not found"
  elif [[ ! -f "$R2C_GZ" ]]; then
    echo "[4/5] skipping r2c PASS filter: $R2C_GZ not found (old svaba run?)"
  else
    echo "[4/5] filter r2c to PASS contigs"
    pass_list="$(mktemp -t "svaba.${ID}.pass_cnames.XXXXXX")"
    # Column positions as documented above: 30=cname, 32=confidence.
    zcat "$BPS_GZ" \
      | awk -F'\t' 'NR>1 && $32=="PASS" {print $30}' \
      | LC_ALL=C sort -u > "$pass_list"
    n_pass=$(wc -l < "$pass_list" | awk '{print $1}')

    if [[ "$n_pass" -eq 0 ]]; then
      echo "      no PASS breakpoints — skipping r2c.pass generation."
      rm -f "$pass_list"
    else
      zcat "$R2C_GZ" \
        | awk -F'\t' -v listf="$pass_list" '
            BEGIN { while ((getline line < listf) > 0) keep[line]=1; close(listf) }
            NR==1                              { print; next }            # header
            $1 == "contig" && ($2 in keep)     { print; next }
            $1 == "read"   && ($2 in keep)     { print; next }
          ' \
        | gzip -c > "$R2C_OUT"
      n_rows=$(zcat "$R2C_OUT" | wc -l | awk '{print $1}')
      echo "      wrote $R2C_OUT  ($n_pass PASS contigs, $n_rows rows incl. header)"
      rm -f "$pass_list"
    fi
  fi
else
  echo "[4/5] skipping r2c PASS filter (--skip-r2c)"
fi

# ========================================================================
# STEP 5 (optional): Split deduped BAMs by QNAME source prefix.
#
# Enabled by --split-by-source (or SPLIT_BY_SOURCE=1). Routes each record
# into ${ID}.${suffix}.${prefix}.bam based on the first 4 chars of QNAME
# (the "source tag" like t001/n001), then sorts+indexes each output.
# Runs for discordant/weird/corrected only (not contigs).
# ========================================================================
if [[ "$SPLIT_BY_SOURCE" == "1" ]]; then
  echo "[5/5] splitting BAMs by QNAME source prefix"
  for suffix in discordant weird corrected; do
    bam="${INPUT_DIR}/${ID}.${suffix}.bam"
    if [[ ! -f "$bam" ]]; then
      echo "      skip: $bam not found"
      continue
    fi
    echo "      splitting $bam"

    TMP_SPLIT=$(mktemp -d "${TMPDIR:-/tmp}/svaba_pp_split.XXXXXXXX")
    header="${TMP_SPLIT}/header.sam"
    "$SAM" view -H "$bam" > "$header"
    # Single-pass route by first 4 chars of QNAME into per-prefix SAM files.
    "$SAM" view "$bam" | awk -v tmpdir="$TMP_SPLIT" '{
      prefix = substr($1, 1, 4)
      out = tmpdir "/body." prefix ".sam"
      print >> out
    }'

    shopt -s nullglob
    bodies=("${TMP_SPLIT}"/body.*.sam)
    shopt -u nullglob
    n_bodies=${#bodies[@]}
    if [[ $n_bodies -gt 0 ]]; then
      per_job_threads=$(( THREADS / n_bodies ))
      [[ $per_job_threads -lt 1 ]] && per_job_threads=1
      max_parallel=$(( THREADS < n_bodies ? THREADS : n_bodies ))
      echo "        reassembling $n_bodies per-prefix BAMs (max_parallel=$max_parallel, threads/job=$per_job_threads)"

      pids=()
      for body in "${bodies[@]}"; do
        while [[ ${#pids[@]} -ge $max_parallel ]]; do
          new_pids=()
          for pid in "${pids[@]}"; do
            kill -0 "$pid" 2>/dev/null && new_pids+=("$pid")
          done
          pids=("${new_pids[@]}")
          [[ ${#pids[@]} -ge $max_parallel ]] && sleep 0.1
        done
        base=$(basename "$body" .sam)
        prefix="${base#body.}"
        out="${OUTPUT_DIR}/${ID}.${suffix}.${prefix}.bam"
        echo "          writing $out"
        (
          cat "$header" "$body" \
            | "$SAM" view -@ "$per_job_threads" -u -bS - \
            | "$SAM" sort -@ "$per_job_threads" -m "$MEM" -o "$out" - \
            && "$SAM" index -@ "$per_job_threads" "$out"
        ) &
        pids+=("$!")
      done
      wait
    fi

    if [[ "${KEEP_TMP:-0}" != "1" ]]; then
      rm -rf "$TMP_SPLIT"
      TMP_SPLIT=""
    fi
  done
else
  echo "[5/5] skipping split-by-source (pass --split-by-source to enable)"
fi

echo "svaba_postprocess.sh: done."
