#!/usr/bin/env bash
#
# sort_and_deduplicate_bps.sh — sort, deduplicate, and PASS-filter a
# svaba bps file.
#
# Reads:   ${ID}.bps.txt.gz   (preferred)
#       or ${ID}.bps.txt      (fallback, uncompressed)
#
# Writes THREE outputs, each sorted by the same key:
#   1. ${ID}.bps.sorted.txt.gz
#        - full sorted list, every row kept.
#   2. ${ID}.bps.sorted.dedup.txt.gz
#        - one row per unique breakpoint pair
#          (chr1,pos1,strand1,chr2,pos2,strand2). When several rows
#          share a breakpoint pair, the one with the HIGHEST maxlod
#          is kept (that's the "best" SV at that junction).
#   3. ${ID}.bps.sorted.dedup.pass.txt.gz
#        - the dedup output further restricted to rows whose
#          confidence column (col 32) equals "PASS".
#
# Preserves the leading "#chr1 ..." header row in every output.
#
# Sort keys, in order:
#   1) chr1    ascending   (natural / version order: chr1, chr2, ..., chrX)
#   2) pos1    ascending
#   3) strand1 ascending
#   4) chr2    ascending
#   5) pos2    ascending
#   6) strand2 ascending
#   7) maxlod  DESCENDING  (higher first; "NA" sorts last)
#
# maxlod is used instead of somlod because the "best" variant at a
# junction is the one most confidently called as a variant at all
# (maxlod = LL(alt) - LL(error)), independent of tumor/normal
# partitioning. somlod can stay low even when maxlod is high, and
# vice versa.
#
# The sort is external — GNU sort spills sorted chunks to a scratch
# directory and merges them, so memory stays bounded even on 10+ GB
# inputs. Works on Linux and macOS; on macOS you'll want GNU coreutils
# installed (`brew install coreutils`) so `gsort` is available.
#
# Usage:
#   ./sort_and_deduplicate_bps.sh <ID>
#   ./sort_and_deduplicate_bps.sh -h | --help
#
# Options via environment variables (all optional):
#   BUFFER_SIZE   per-chunk memory for sort       (default: 2G)
#   PARALLEL      threads used by the sort phase  (default: 4)
#   INPUT_DIR     where to look for the input     (default: .)
#   OUTPUT_DIR    where to write the outputs      (default: $INPUT_DIR)
#   KEEP_TMP      1 = keep the scratch dir for debugging (default: unset)
#
# Examples:
#   sort_and_deduplicate_bps.sh tumor_vs_normal
#   BUFFER_SIZE=8G PARALLEL=8 sort_and_deduplicate_bps.sh tumor_vs_normal
#   INPUT_DIR=/data/run42 OUTPUT_DIR=/tmp sort_and_deduplicate_bps.sh run42
#

set -euo pipefail

# ------------------------------------------------------------------ usage ---

print_usage() {
    cat <<'EOF'
sort_and_deduplicate_bps.sh — sort + dedup + PASS-filter a svaba bps file

Usage:
  sort_and_deduplicate_bps.sh <ID>
  sort_and_deduplicate_bps.sh -h | --help

Reads  <INPUT_DIR>/<ID>.bps.txt.gz  (falls back to <ID>.bps.txt if no .gz)
Writes (all under <OUTPUT_DIR>):
  <ID>.bps.sorted.txt.gz             — full sorted list
  <ID>.bps.sorted.dedup.txt.gz       — best SV per breakpoint pair
  <ID>.bps.sorted.dedup.pass.txt.gz  — dedup output restricted to PASS

Sort keys (in order):
  1. chr1     ascending   (natural / version order)
  2. pos1     ascending
  3. strand1  ascending
  4. chr2     ascending
  5. pos2     ascending
  6. strand2  ascending
  7. maxlod   DESCENDING  (higher first; "NA" sorts last)

Dedup:
  One row kept per unique (chr1,pos1,strand1,chr2,pos2,strand2) tuple.
  Because we sort maxlod DESCENDING before deduping, the row kept is the
  one with the highest maxlod at that junction — the "best" SV.

PASS filter:
  Keeps only rows whose confidence column (col 32) equals literal "PASS".

Options (environment variables):
  BUFFER_SIZE   per-chunk memory for sort       (default: 2G)
  PARALLEL      threads used by the sort phase  (default: 4)
  INPUT_DIR     where to look for the input     (default: .)
  OUTPUT_DIR    where to write the outputs      (default: $INPUT_DIR)
  KEEP_TMP      1 = keep the scratch dir for debugging (default: unset)

Notes:
  * Uses GNU sort. On macOS install it via `brew install coreutils` —
    the script will pick up `gsort` automatically. BSD sort does not
    support -V or --parallel.
  * Memory is bounded by BUFFER_SIZE; scratch disk usage is roughly the
    uncompressed size of the input (often 6-10x the .gz size).
EOF
}

# ---------------------------------------------------------------- args -----

if [[ $# -eq 0 ]]; then
    print_usage
    exit 1
fi

case "${1:-}" in
    -h|--help)
        print_usage
        exit 0
        ;;
    -*)
        echo "sort_and_deduplicate_bps.sh: unknown option '$1'" >&2
        echo "Run 'sort_and_deduplicate_bps.sh --help' for usage." >&2
        exit 2
        ;;
esac

ID=$1
if [[ -z "$ID" ]]; then
    echo "sort_and_deduplicate_bps.sh: <ID> is required" >&2
    print_usage
    exit 2
fi

INPUT_DIR=${INPUT_DIR:-.}
OUTPUT_DIR=${OUTPUT_DIR:-$INPUT_DIR}
BUFFER_SIZE=${BUFFER_SIZE:-2G}
PARALLEL=${PARALLEL:-4}

# ----------------------------------------------------- input detection -----
#
# Prefer the .gz. Fall back to plain .txt. Bail if neither is readable.

IN_GZ="${INPUT_DIR%/}/${ID}.bps.txt.gz"
IN_TXT="${INPUT_DIR%/}/${ID}.bps.txt"

if   [[ -r "$IN_GZ"  ]]; then IN="$IN_GZ";  READER="gzip -dc"
elif [[ -r "$IN_TXT" ]]; then IN="$IN_TXT"; READER="cat"
else
    echo "sort_and_deduplicate_bps.sh: cannot read input — tried:" >&2
    echo "  $IN_GZ"  >&2
    echo "  $IN_TXT" >&2
    exit 3
fi

OUT_SORTED="${OUTPUT_DIR%/}/${ID}.bps.sorted.txt.gz"
OUT_DEDUP="${OUTPUT_DIR%/}/${ID}.bps.sorted.dedup.txt.gz"
OUT_PASS="${OUTPUT_DIR%/}/${ID}.bps.sorted.dedup.pass.txt.gz"

mkdir -p "$(dirname "$OUT_SORTED")"

# ----------------------------------------------------------- tool lookup ---

# Prefer `gsort` (GNU coreutils on macOS via homebrew); fall back to `sort`
# if it reports itself as GNU.
pick_sort() {
    if command -v gsort >/dev/null 2>&1; then
        echo gsort
        return 0
    fi
    if sort --version 2>/dev/null | head -n 1 | grep -qi 'gnu'; then
        echo sort
        return 0
    fi
    return 1
}

if ! SORT_BIN=$(pick_sort); then
    echo "sort_and_deduplicate_bps.sh: GNU sort not found." >&2
    echo "  On macOS: brew install coreutils   (provides 'gsort')" >&2
    echo "  On Linux: install the 'coreutils' package" >&2
    exit 4
fi

# ---------------------------------------------------------- scratch dir ---

TMP=$(mktemp -d "${TMPDIR:-/tmp}/svaba_bps_sortdedup.XXXXXXXX")

cleanup() {
    if [[ "${KEEP_TMP:-0}" == "1" ]]; then
        echo "sort_and_deduplicate_bps.sh: leaving scratch dir at $TMP" >&2
    else
        rm -rf "$TMP"
    fi
}
trap cleanup EXIT INT TERM

# ---------------------------------------------- columns in the bps file ---
#
# Hard-coded here so awk stays simple. If the svaba output format ever
# changes column positions, update these two numbers.
#
#   col 32 = confidence    (e.g. PASS, LOWMAPQ, BLACKLIST, LOWICSUPPORT)
#   col 38 = maxlod        (numeric, can be "NA")
#
# Sort keys 1..6 below correspond to chr1,pos1,strand1,chr2,pos2,strand2
# at cols 1,2,3,4,5,6.

CONF_COL=32
MAXLOD_COL=38

# ========================================================================
# STAGE 1: sort
# ========================================================================

echo "[1/3] sorting  $IN -> $OUT_SORTED" >&2
echo "      sort=$SORT_BIN buffer=$BUFFER_SIZE parallel=$PARALLEL tmp=$TMP" >&2

# Peek at the first line to split the header from the body without two
# passes over the input.
FIRST_LINE=$($READER "$IN" | head -n 1 || true)
HAS_HEADER=0
if [[ "${FIRST_LINE:0:1}" == "#" ]]; then
    HAS_HEADER=1
fi

{
    if [[ $HAS_HEADER -eq 1 ]]; then
        printf '%s\n' "$FIRST_LINE"
        $READER "$IN" | tail -n +2
    else
        $READER "$IN"
    fi | {
        if [[ $HAS_HEADER -eq 1 ]]; then
            # Pass the header through untouched; sort only the body.
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
            -S "$BUFFER_SIZE" \
            --parallel="$PARALLEL" \
            -T "$TMP"
    }
} | gzip -c > "$OUT_SORTED"

# ========================================================================
# STAGE 2: dedup by breakpoint pair (keep first = highest maxlod)
# ========================================================================
#
# Because stage 1 already sorted by (chr1,pos1,strand1,chr2,pos2,strand2,
# -maxlod), the first row awk sees for each breakpoint-pair key is the
# row with the highest maxlod at that junction. We just print the first
# occurrence of each key.

echo "[2/3] deduping $OUT_SORTED -> $OUT_DEDUP" >&2

gzip -dc "$OUT_SORTED" | awk -F'\t' -v OFS='\t' '
    # header: print untouched, do not contribute to dedup
    NR == 1 && $0 ~ /^#/ { print; next }

    {
        key = $1 FS $2 FS $3 FS $4 FS $5 FS $6
        if (!(key in seen)) {
            seen[key] = 1
            print
        }
    }
' | gzip -c > "$OUT_DEDUP"

# ========================================================================
# STAGE 3: PASS-only filter
# ========================================================================

echo "[3/3] PASS-filtering $OUT_DEDUP -> $OUT_PASS" >&2

gzip -dc "$OUT_DEDUP" | awk -F'\t' -v OFS='\t' -v CCOL="$CONF_COL" '
    # header: print untouched
    NR == 1 && $0 ~ /^#/ { print; next }

    # keep only PASS rows
    $CCOL == "PASS" { print }
' | gzip -c > "$OUT_PASS"

# ------------------------------------------------------- final summary ---

n_lines() {
    # counts data rows (excludes the single "#" header row if present)
    local f=$1
    local total header
    total=$(gzip -dc "$f" | wc -l)
    header=$(gzip -dc "$f" | head -n 1 | grep -c '^#' || true)
    echo $(( total - header ))
}

echo "sort_and_deduplicate_bps.sh: done." >&2
echo "  sorted             : $OUT_SORTED  ($(n_lines "$OUT_SORTED") rows)" >&2
echo "  sorted + dedup     : $OUT_DEDUP   ($(n_lines "$OUT_DEDUP") rows)" >&2
echo "  sorted + dedup+PASS: $OUT_PASS    ($(n_lines "$OUT_PASS") rows)" >&2
