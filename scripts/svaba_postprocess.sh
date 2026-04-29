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
#        ${ID}.bps.sorted.dedup.pass.somatic.txt.gz — PASS rows that are also
#                                              somatic (somlod >= 0 for v3
#                                              format, somatic_state == 1 for
#                                              legacy format)
#
#   4. Filters ${ID}.r2c.txt.gz down to PASS breakpoints — writes both
#      ${ID}.r2c.pass.txt.gz (all PASS calls) and
#      ${ID}.r2c.pass.somatic.txt.gz (PASS calls with somlod >= 1,
#      i.e. the somatic subset) for interactive use in bps_explorer.html
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
#       --skip-dedup      keep step 2 but tell svaba postprocess to skip the
#                         dedup substep (passes --sort-only through). Useful
#                         for quick reruns of sort + @PG stamp + index
#                         without redoing the expensive dedup pass. Sort is
#                         also auto-skipped when the BAM's header already
#                         declares @HD SO:coordinate, so a repeat rerun is
#                         effectively a no-op except for index rebuild.
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
SKIP_DEDUP=0
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
    --skip-dedup)      SKIP_DEDUP=1; shift ;;
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
echo "                    : split_by_source=$SPLIT_BY_SOURCE skip_bam=$SKIP_BAM skip_dedup=$SKIP_DEDUP skip_bps=$SKIP_BPS skip_r2c=$SKIP_R2C"

# Useful for absolute-path filenames in the steps below. Use trailing-slash-
# free form so we always concatenate `${DIR}/${name}`.
INPUT_DIR=${INPUT_DIR%/}
OUTPUT_DIR=${OUTPUT_DIR%/}

# ========================================================================
# STEP 1: Merge per-thread BAMs and per-thread r2c.txt.gz.
#
# svaba run with multiple threads produces ${ID}.thread*.${suffix}.bam for
# suffix in {discordant,weird,corrected}, and (with --dump-reads)
# ${ID}.thread*.r2c.txt.gz. This step coalesces each family into a single
# ${ID}.${suffix}.bam / ${ID}.r2c.txt.gz.
#
# Per-thread BAMs are merged with `samtools merge` — standard BAM merge
# with header reconciliation.
#
# Per-thread r2c.txt.gz files are merged with plain `cat`: gzip is a
# concatenation-safe format (RFC 1952), so concatenating multiple .gz
# streams produces a valid .gz whose decompressed output is the
# concatenation of each member's decompressed output. svaba guarantees
# only the first worker (threadId == 1; the worker pool numbers
# threads 1..N) emits the TSV column-header line, and the awk|sort
# below places that file first in the merge, so the merged file has
# exactly one header at the top regardless of thread count.
#
# Single-file inputs (only one thread) are renamed. No-file inputs are
# silently skipped — this happens when the svaba run was single-threaded,
# or when --dump-reads wasn't set (no r2c files to merge).
# ========================================================================
if [[ $SKIP_BAM -eq 0 ]]; then
  echo "[1/5] merging per-thread BAMs"
  for suffix in discordant weird corrected; do
    target="${INPUT_DIR}/${ID}.${suffix}.bam"
    shopt -s nullglob
    bam_files=("${INPUT_DIR}/${ID}".thread*."${suffix}".bam)
    shopt -u nullglob
    if [[ ${#bam_files[@]} -gt 1 ]]; then
      echo "      merge ${#bam_files[@]} thread BAMs for '$suffix' -> $target"
      "$SAM" merge -f -@ "$THREADS" "$target" "${bam_files[@]}"
      rm -f "${bam_files[@]}"
    elif [[ ${#bam_files[@]} -eq 1 ]]; then
      echo "      mv ${bam_files[0]} -> $target"
      mv "${bam_files[0]}" "$target"
    fi
  done

  # Merge per-thread r2c.txt.gz → ${ID}.r2c.txt.gz. gzip is concat-safe,
  # so `cat a.gz b.gz > c.gz` is valid. Thread 1's file is first in the
  # numeric sort order (worker pool numbers threads 1..N; see
  # threadpool.h), and thread 1 is the one that emits the TSV header,
  # so its header line ends up at the top of the merged file.
  target="${INPUT_DIR}/${ID}.r2c.txt.gz"
  shopt -s nullglob
  r2c_files=("${INPUT_DIR}/${ID}".thread*.r2c.txt.gz)
  shopt -u nullglob
  if [[ ${#r2c_files[@]} -gt 1 ]]; then
    echo "      merge ${#r2c_files[@]} thread r2c files -> $target"
    # Sort by thread index so thread 1 (with the header) is first.
    # Filenames are ${ID}.threadN.r2c.txt.gz — awk pulls the numeric
    # thread index out, sort -k1,1n orders it numerically, cut drops
    # the index so we get just the filenames back in sorted order.
    #
    # NB: avoid `mapfile` / `readarray` — they're bash 4+, and macOS
    # ships bash 3.2 as /bin/bash, which `#!/usr/bin/env bash` often
    # resolves to. The `while read` loop below is fully portable.
    r2c_sorted=()
    while IFS= read -r f; do
      r2c_sorted+=("$f")
    done < <(
      printf '%s\n' "${r2c_files[@]}" \
      | awk -F'.thread|.r2c.txt.gz' '{print $(NF-1)"\t"$0}' \
      | sort -k1,1n \
      | cut -f2-
    )
    cat "${r2c_sorted[@]}" > "$target"
    rm -f "${r2c_sorted[@]}"
  elif [[ ${#r2c_files[@]} -eq 1 ]]; then
    echo "      mv ${r2c_files[0]} -> $target"
    mv "${r2c_files[0]}" "$target"
  fi
else
  echo "[1/5] skipping BAM + r2c merge (--skip-bam)"
fi

# ========================================================================
# STEP 2: svaba postprocess — sort, stream-dedup, stamp @PG, build .bai.
#
# Heavy lifting lives in the C++ subcommand. We just pass it the thread /
# memory budget; svaba postprocess splits THREADS across suffixes and passes
# -m to samtools sort per-thread.
#
# --skip-dedup at this layer maps to `--sort-only` on the C++ CLI, which
# still runs the sort + @PG-stamp + index phases but skips the
# stream-dedup pass. Combined with the C++'s auto-skip-sort-when-already-
# sorted check (via @HD SO:coordinate), a rerun with --skip-dedup on
# already-postprocessed files is effectively a no-op other than
# re-stamping PG + rebuilding the .bai.
# ========================================================================
if [[ $SKIP_BAM -eq 0 ]]; then
  pp_args=(-i "$ID" -t "$THREADS" -m "$MEM")
  if [[ $SKIP_DEDUP -eq 1 ]]; then
    pp_args+=(--sort-only)
    dedup_banner="(dedup skipped via --skip-dedup)"
  else
    dedup_banner=""
  fi
  echo "[2/5] svaba postprocess: sort + dedup + @PG + BAI index (threads=$THREADS, mem=$MEM) $dedup_banner"
  # Indent all sub-output and suppress the harmless htslib "idx_find_and_load"
  # warnings (index files don't exist yet when the reader opens the BAM —
  # they're built as the final step of each suffix's pipeline). Also strip
  # the NOTE line about those warnings since we're filtering them here.
  # Use awk instead of grep -v so a zero-match case doesn't trigger pipefail.
  (cd "$INPUT_DIR" && "$SVABA" postprocess "${pp_args[@]}" 2>&1) \
    | awk '!/idx_find_and_load/ && !/Could not retrieve index file.*expected and can be ignored/' \
    | sed 's/^/      /'
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
OUT_PASS_SOM="${OUTPUT_DIR}/${ID}.bps.sorted.dedup.pass.somatic.txt.gz"

# Column positions — auto-detected from the header when available (see
# detect_columns below). These are the v3 defaults; the old format has
# different positions (e.g. confidence at col 26 instead of 32).
CONF_COL=32
SOMATIC_COL=36
SOMLOD_COL=37
MAXLOD_COL=38
CNAME_COL=30

pick_sort() {
  if command -v gsort >/dev/null 2>&1; then
    echo gsort; return 0
  fi
  if sort --version 2>/dev/null | head -n 1 | grep -qi 'gnu'; then
    echo sort; return 0
  fi
  return 1
}

# Auto-detect column positions from the header line. Looks for known
# column names and sets CONF_COL, SOMATIC_COL, SOMLOD_COL, MAXLOD_COL.
# Falls back to v3 defaults if the header is missing or unrecognised.
detect_columns() {
  local header_line="$1"
  # strip leading # if present
  local clean="${header_line#\#}"

  # Use awk to find column indices by name. Tries both old and new
  # column names (e.g. "somatic_score" vs the unnamed v3 col,
  # "contig" vs "cname").
  eval "$(echo "$clean" | awk -F'\t' '{
    for (i = 1; i <= NF; i++) {
      col = $i
      # strip whitespace
      gsub(/^[ \t]+|[ \t]+$/, "", col)
      # normalise to lowercase for matching
      lcol = tolower(col)

      if (lcol == "confidence")          printf "CONF_COL=%d\n", i
      if (lcol == "somatic_score")       printf "SOMATIC_COL=%d\n", i
      if (lcol == "somatic_lod")         printf "SOMLOD_COL=%d\n", i
      if (lcol == "tlod" || lcol == "maxlod" || lcol == "max_lod")
                                         printf "MAXLOD_COL=%d\n", i
      if (lcol == "contig" || lcol == "cname")
                                         printf "CNAME_COL=%d\n", i
    }
  }')"
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

    # Auto-detect column positions from the header.
    if [[ $HAS_HEADER -eq 1 ]]; then
      detect_columns "$FIRST_LINE"
    fi
    echo "      columns: confidence=$CONF_COL cname=$CNAME_COL somatic=$SOMATIC_COL somlod=$SOMLOD_COL maxlod=$MAXLOD_COL"

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

    # --- 3d) PASS + somatic filter ---
    # Use somatic_score == "1" as the sole somatic gate. The old format's
    # "somatic_lod" column is actually n.LO_n (normal variant LO) — high
    # values mean strong normal evidence, which argues AGAINST somatic.
    # The v3 "somlod" is a proper somatic log-odds, but somatic_score is
    # set correctly in both formats, so we use it uniformly.
    gzip -dc "$OUT_PASS" | awk -F'\t' -v OFS='\t' \
        -v SCOL="$SOMATIC_COL" '
      NR == 1 && $0 ~ /^#/ { print; next }
      $SCOL == "1" { print }
    ' | gzip -c > "$OUT_PASS_SOM"
    echo "      wrote $OUT_PASS_SOM"
  fi
else
  echo "[3/5] skipping bps sort/dedup/PASS (--skip-bps)"
fi

# ========================================================================
# STEP 4: Filter r2c.txt.gz to PASS contigs (+ PASS-somatic subset).
#
# The structured r2c TSV emitted by `svaba run` contains one row per
# assembled contig + one row per r2c-aligned read, for every variant-
# bearing contig. At full genome scale this can be tens of millions of
# rows — too heavy for the interactive viewer. Here we produce TWO
# filtered outputs from a single pass over r2c.txt.gz:
#
#   ${ID}.r2c.pass.txt.gz          - rows for all PASS variants
#   ${ID}.r2c.pass.somatic.txt.gz  - rows for PASS variants that are also
#                                    somatic (somlod >= 1 for v3 format,
#                                    somatic_score == 1 for legacy)
#
# Both are keyed on cname (col 30 in bps.txt.gz, col 2 in r2c.txt.gz);
# contig and read rows survive together since they share col 2.
#
# The somatic subset is a strict subset of the pass set, and both are
# written by the same awk invocation over r2c.txt.gz via awk's pipe-to-
# command output — one decompress+iterate pass, two gzip compressors
# running in parallel as child processes. Avoids reading the multi-GB
# r2c file twice.
# ========================================================================
R2C_GZ="${INPUT_DIR}/${ID}.r2c.txt.gz"
R2C_OUT="${OUTPUT_DIR}/${ID}.r2c.pass.txt.gz"
R2C_OUT_SOM="${OUTPUT_DIR}/${ID}.r2c.pass.somatic.txt.gz"

if [[ $SKIP_R2C -eq 0 ]]; then
  if [[ ! -f "$BPS_GZ" ]]; then
    echo "[4/5] skipping r2c PASS filter: $BPS_GZ not found"
  elif [[ ! -f "$R2C_GZ" ]]; then
    echo "[4/5] skipping r2c PASS filter: $R2C_GZ not found"
    echo "       (svaba run must be invoked with --dump-reads to produce it;"
    echo "        or this is an older svaba run from before r2c.txt.gz existed.)"
  else
    echo "[4/5] filter r2c to PASS contigs (+ PASS-somatic subset)"
    pass_list="$(mktemp -t "svaba.${ID}.pass_cnames.XXXXXX")"
    som_list="$(mktemp  -t "svaba.${ID}.som_cnames.XXXXXX")"
    # Column positions are auto-detected from the header in step 3 above
    # (CONF_COL, CNAME_COL, SOMLOD_COL, SOMATIC_COL). If step 3 was
    # skipped, re-detect from bps.txt.gz's own header here.
    if [[ $SKIP_BPS -eq 1 ]]; then
      BPS_HEADER=$(gzip -dc "$BPS_GZ" | head -n 1 || true)
      [[ "${BPS_HEADER:0:1}" == "#" ]] && detect_columns "$BPS_HEADER"
      echo "      columns: confidence=$CONF_COL cname=$CNAME_COL somatic=$SOMATIC_COL somlod=$SOMLOD_COL"
    fi
    # One awk pass builds both cname sets so we only decompress bps.txt.gz
    # once. NB: use `gzip -dc` (not `zcat`) so this works on macOS — BSD
    # zcat hunts for a .Z suffix rather than .gz.
    # Somatic criterion: somatic_score == "1" (reliable in both old and v3 formats).
    gzip -dc "$BPS_GZ" \
      | awk -F'\t' -v pass_f="$pass_list" -v som_f="$som_list" \
            -v CCOL="$CONF_COL" -v NCOL="$CNAME_COL" \
            -v SCOL="$SOMATIC_COL" '
          NR>1 && $CCOL=="PASS" {
            print $NCOL >> pass_f
            if ($SCOL == "1")
              { print $NCOL >> som_f }
          }
        '
    LC_ALL=C sort -u "$pass_list" -o "$pass_list"
    LC_ALL=C sort -u "$som_list"  -o "$som_list"
    n_pass=$(wc -l < "$pass_list" | awk '{print $1}')
    n_som=$(wc  -l < "$som_list"  | awk '{print $1}')

    if [[ "$n_pass" -eq 0 ]]; then
      echo "      no PASS breakpoints — skipping both r2c.pass and r2c.pass.somatic."
      rm -f "$pass_list" "$som_list"
    else
      # Single-pass r2c filter writing both outputs via awk's pipe-to-
      # command feature. The two `gzip -c > file` child processes run in
      # parallel, so compression of the two outputs happens concurrently
      # on separate cores. Each output gets the same header row and then
      # whichever filtered rows match its cname set.
      gzip -dc "$R2C_GZ" \
        | awk -F'\t' -v pass_f="$pass_list"   -v som_f="$som_list" \
                     -v pass_out="$R2C_OUT"   -v som_out="$R2C_OUT_SOM" '
            BEGIN {
              while ((getline line < pass_f) > 0) keep_pass[line] = 1
              close(pass_f)
              while ((getline line < som_f)  > 0) keep_som[line]  = 1
              close(som_f)
              pass_cmd = "gzip -c > " pass_out
              som_cmd  = "gzip -c > " som_out
            }
            NR == 1 {
              print | pass_cmd
              print | som_cmd
              next
            }
            ($1 == "contig" || $1 == "read") {
              if ($2 in keep_pass) print | pass_cmd
              if ($2 in keep_som)  print | som_cmd
            }
            END {
              close(pass_cmd)
              close(som_cmd)
            }
          '

      n_pass_rows=$(gzip -dc "$R2C_OUT"     | wc -l | awk '{print $1}')
      n_som_rows=$(gzip  -dc "$R2C_OUT_SOM" | wc -l | awk '{print $1}')
      echo "      wrote $R2C_OUT        ($n_pass PASS contigs, $n_pass_rows rows incl. header)"
      echo "      wrote $R2C_OUT_SOM  ($n_som PASS+somatic [somlod>=1] contigs, $n_som_rows rows incl. header)"
      rm -f "$pass_list" "$som_list"
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
