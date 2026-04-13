#!/usr/bin/env bash
#
# sort_bps.sh — sort a svaba bps.txt.gz file by (chr1, pos1, chr2, pos2, -somlod)
#
# Reads:   ${ID}.bps.txt.gz
# Writes:  ${ID}.bps.sorted.txt.gz
#
# Preserves the leading "#chr1 ..." header row, then sorts the body by:
#   1) chr1 (natural/version order)   ascending
#   2) pos1 (numeric)                 ascending
#   3) chr2 (natural/version order)   ascending
#   4) pos2 (numeric)                 ascending
#   5) somlod (general numeric)       DESCENDING (higher first)
#
# The sort is external — GNU sort spills sorted chunks to a scratch directory
# and merges them, so memory stays bounded even on 10+ GB inputs. The scratch
# directory is auto-created with mktemp and cleaned up on exit. Works on both
# Linux and macOS; on macOS you'll want GNU coreutils installed (`brew install
# coreutils`) so `gsort` is available — BSD sort does not support -V / --parallel.
#
# Usage:
#   ./sort_bps.sh <ID>
#   ./sort_bps.sh -h | --help
#
# Options via environment variables (all optional):
#   BUFFER_SIZE   per-chunk memory for sort       (default: 2G)
#   PARALLEL      threads used by the sort phase  (default: 4)
#   INPUT_DIR     where to look for the input     (default: .)
#   OUTPUT_DIR    where to write the output       (default: $INPUT_DIR)
#   KEEP_TMP      1 = don't delete the scratch dir on exit (for debugging)
#
# Examples:
#   ./sort_bps.sh tumor_vs_normal
#   BUFFER_SIZE=8G PARALLEL=8 ./sort_bps.sh tumor_vs_normal
#   INPUT_DIR=/data/run42 OUTPUT_DIR=/tmp ./sort_bps.sh run42
#

set -euo pipefail

# ------------------------------------------------------------------ usage ---

print_usage() {
    cat <<'EOF'
sort_bps.sh — sort a svaba bps.txt.gz file

Usage:
  sort_bps.sh <ID>
  sort_bps.sh -h | --help

Reads  <INPUT_DIR>/<ID>.bps.txt.gz
Writes <OUTPUT_DIR>/<ID>.bps.sorted.txt.gz

Sort keys (in order):
  1. chr1  ascending   (natural / version order — chr1, chr2, ..., chrX)
  2. pos1  ascending
  3. chr2  ascending
  4. pos2  ascending
  5. somlod DESCENDING (higher somatic LOD first; "NA" sorts last)

Options (environment variables):
  BUFFER_SIZE   per-chunk memory for sort       (default: 2G)
  PARALLEL      threads used by the sort phase  (default: 4)
  INPUT_DIR     where to look for the input     (default: .)
  OUTPUT_DIR    where to write the output       (default: $INPUT_DIR)
  KEEP_TMP      1 = keep the scratch dir for debugging (default: unset)

Examples:
  sort_bps.sh tumor_vs_normal
  BUFFER_SIZE=8G PARALLEL=8 sort_bps.sh tumor_vs_normal
  INPUT_DIR=/data/run42 OUTPUT_DIR=/tmp sort_bps.sh run42

Notes:
  * Uses GNU sort. On macOS install it via `brew install coreutils` — the
    script will pick up `gsort` automatically. BSD sort does not support -V
    or --parallel and will not work.
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
        echo "sort_bps.sh: unknown option '$1'" >&2
        echo "Run 'sort_bps.sh --help' for usage." >&2
        exit 2
        ;;
esac

ID=$1
if [[ -z "$ID" ]]; then
    echo "sort_bps.sh: <ID> is required" >&2
    print_usage
    exit 2
fi

INPUT_DIR=${INPUT_DIR:-.}
OUTPUT_DIR=${OUTPUT_DIR:-$INPUT_DIR}
BUFFER_SIZE=${BUFFER_SIZE:-2G}
PARALLEL=${PARALLEL:-4}

IN="${INPUT_DIR%/}/${ID}.bps.txt.gz"
OUT="${OUTPUT_DIR%/}/${ID}.bps.sorted.txt.gz"

if [[ ! -r "$IN" ]]; then
    echo "sort_bps.sh: cannot read input '$IN'" >&2
    exit 3
fi

mkdir -p "$(dirname "$OUT")"

# ----------------------------------------------------------- tool lookup ---

# Prefer `gsort` (GNU coreutils on macOS via homebrew); fall back to `sort`
# if it reports itself as GNU. Bail out if neither works.
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
    echo "sort_bps.sh: GNU sort not found." >&2
    echo "  On macOS: brew install coreutils   (provides 'gsort')" >&2
    echo "  On Linux: install the 'coreutils' package" >&2
    exit 4
fi

# zcat on macOS appends .Z instead of reading .gz — use `gzip -dc` to stay
# portable.
ZCAT="gzip -dc"

# ---------------------------------------------------------- scratch dir ---

# mktemp -d is portable across Linux and macOS. The template form with a
# trailing XXXXXX works on both.
TMP=$(mktemp -d "${TMPDIR:-/tmp}/svaba_bps_sort.XXXXXXXX")

cleanup() {
    if [[ "${KEEP_TMP:-0}" == "1" ]]; then
        echo "sort_bps.sh: leaving scratch dir at $TMP" >&2
    else
        rm -rf "$TMP"
    fi
}
trap cleanup EXIT INT TERM

# ------------------------------------------------------------- the sort ---

echo "sort_bps.sh: sorting $IN -> $OUT" >&2
echo "sort_bps.sh: sort=$SORT_BIN buffer=$BUFFER_SIZE parallel=$PARALLEL tmp=$TMP" >&2

# Extract the header once (if present) so we don't have to zcat twice.
HEADER_FILE="$TMP/header.txt"
BODY_PIPE=0

# Peek at the first line. If it starts with '#', save it; otherwise there is
# no header and we just sort the whole file.
FIRST_LINE=$($ZCAT "$IN" | head -n 1 || true)
if [[ "${FIRST_LINE:0:1}" == "#" ]]; then
    printf '%s\n' "$FIRST_LINE" > "$HEADER_FILE"
    BODY_PIPE=1
fi

{
    if [[ $BODY_PIPE -eq 1 ]]; then
        cat "$HEADER_FILE"
        $ZCAT "$IN" | tail -n +2
    else
        $ZCAT "$IN"
    fi | {
        # The header (if any) is the first line; pass it through untouched,
        # then pipe the rest into sort. Using `read` here keeps everything
        # streaming — no intermediate uncompressed file on disk.
        if [[ $BODY_PIPE -eq 1 ]]; then
            IFS= read -r hdr
            printf '%s\n' "$hdr"
        fi
        LC_ALL=C "$SORT_BIN" \
            -t $'\t' \
            -k1,1V -k2,2n -k4,4V -k5,5n -k37,37gr \
            -S "$BUFFER_SIZE" \
            --parallel="$PARALLEL" \
            -T "$TMP"
    }
} | gzip -c > "$OUT"

echo "sort_bps.sh: done -> $OUT" >&2
