#!/usr/bin/env bash
#
# sort_and_dedupe_bps_old.sh — sort + dedup + PASS-filter for legacy
# (pre-2.0) svaba bps.txt.gz files.
#
# Same logic as the v3 postprocess step but with auto-detected column
# positions from the header row. Works on any bps.txt.gz schema as
# long as the header has #chr1, confidence, and either maxlod or tlod.
#
# Usage:
#   sort_and_dedupe_bps_old.sh [options] <input.bps.txt.gz>
#
#   -o DIR         output directory (default: same as input)
#   -t THREADS     parallel sort threads (default: 4)
#   -S SIZE        sort buffer size (default: 4G)
#   -h, --help     this message
#
# Outputs:
#   <basename>.sorted.txt.gz
#   <basename>.sorted.dedup.txt.gz
#   <basename>.sorted.dedup.pass.txt.gz
#

set -exo pipefail

THREADS=4
SORT_BUFFER="4G"
OUT_DIR=""

print_usage() { sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; }

while [[ $# -gt 0 ]]; do
  case "$1" in
    -o)        OUT_DIR="$2";      shift 2 ;;
    -t)        THREADS="$2";      shift 2 ;;
    -S)        SORT_BUFFER="$2";  shift 2 ;;
    -h|--help) print_usage;       exit 0 ;;
    -*)        echo "unknown option: $1" >&2; exit 2 ;;
    *)         break ;;
  esac
done

if [[ $# -lt 1 ]]; then
  echo "usage: sort_and_dedupe_bps_old.sh [options] <input.bps.txt.gz>" >&2
  exit 2
fi

INPUT="$1"
if [[ ! -r "$INPUT" ]]; then
  echo "cannot read: $INPUT" >&2
  exit 1
fi

# Determine reader
case "$INPUT" in
  *.gz)  READER="gzip -dc" ;;
  *)     READER="cat" ;;
esac

# Derive output names
BASENAME=$(basename "$INPUT" .txt.gz)
BASENAME=$(basename "$BASENAME" .txt)
[[ -z "$OUT_DIR" ]] && OUT_DIR=$(dirname "$INPUT")

OUT_SORTED="${OUT_DIR}/${BASENAME}.sorted.txt.gz"
OUT_DEDUP="${OUT_DIR}/${BASENAME}.sorted.dedup.txt.gz"
OUT_PASS="${OUT_DIR}/${BASENAME}.sorted.dedup.pass.txt.gz"

# Pick GNU sort
if command -v gsort >/dev/null 2>&1; then
  SORT_BIN=gsort
elif sort --version 2>&1 | grep -qi 'gnu' 2>/dev/null; then
  SORT_BIN=sort
else
  echo "GNU sort not found (macOS: 'brew install coreutils' for gsort)" >&2
  exit 1
fi

# --- Auto-detect column positions from header ---
HEADER=$($READER "$INPUT" | head -n 1 || true)
if [[ "${HEADER:0:1}" != "#" ]]; then
  echo "ERROR: first line doesn't start with '#' — not a bps.txt.gz header?" >&2
  echo "  got: ${HEADER:0:80}..." >&2
  exit 1
fi

# Split header into array, find columns by name
find_col() {
  local name="$1"
  local col
  col=$(echo "$HEADER" | awk -F'\t' -v name="$name" '{
    for (i=1; i<=NF; i++) if ($i == name) { print i; exit }
  }')
  if [[ -z "$col" ]]; then
    return 1
  fi
  echo "$col"
}

CONF_COL=$(find_col "confidence") || { echo "ERROR: 'confidence' column not found in header" >&2; exit 1; }

# Try maxlod first (v3), then tlod (legacy)
if LOD_COL=$(find_col "maxlod" 2>/dev/null); then
  LOD_NAME="maxlod"
elif LOD_COL=$(find_col "tlod" 2>/dev/null); then
  LOD_NAME="tlod"
else
  echo "ERROR: neither 'maxlod' nor 'tlod' column found in header" >&2
  exit 1
fi

echo "=== sort_and_dedupe_bps_old.sh ==="
echo "  input:      $INPUT"
echo "  confidence: col $CONF_COL"
echo "  ${LOD_NAME}:     col $LOD_COL"
echo "  sort:       $SORT_BIN  buffer=$SORT_BUFFER  parallel=$THREADS"
echo ""

TMP_DIR=$(mktemp -d "${TMPDIR:-/tmp}/svaba_sort_bps.XXXXXXXX")
trap 'rm -rf "$TMP_DIR"' EXIT INT TERM

# --- 1) sort ---
echo "[1/3] sorting..."
{
  # Emit header first, then sort the data rows
  printf '%s\n' "$HEADER"
  $READER "$INPUT" | tail -n +2 | \
    LC_ALL=C "$SORT_BIN" \
      -t $'\t' \
      -k1,1V \
      -k2,2n \
      -k3,3 \
      -k4,4V \
      -k5,5n \
      -k6,6 \
      -k${LOD_COL},${LOD_COL}gr \
      -S "$SORT_BUFFER" \
      --parallel="$THREADS" \
      -T "$TMP_DIR"
} | gzip -c > "$OUT_SORTED"
echo "  wrote $OUT_SORTED"

# --- 2) dedup by breakpoint pair (first row per key wins = highest lod) ---
echo "[2/3] deduplicating..."
gzip -dc "$OUT_SORTED" | awk -F'\t' -v OFS='\t' '
  NR == 1 && $0 ~ /^#/ { print; next }
  { key = $1 FS $2 FS $3 FS $4 FS $5 FS $6
    if (!(key in seen)) { seen[key] = 1; print } }
' | gzip -c > "$OUT_DEDUP"
echo "  wrote $OUT_DEDUP"

# --- 3) PASS-filter ---
echo "[3/3] PASS-filtering..."
gzip -dc "$OUT_DEDUP" | awk -F'\t' -v OFS='\t' -v CCOL="$CONF_COL" '
  NR == 1 && $0 ~ /^#/ { print; next }
  $CCOL == "PASS" { print }
' | gzip -c > "$OUT_PASS"
echo "  wrote $OUT_PASS"

# Summary
N_SORTED=$(gzip -dc "$OUT_SORTED" | tail -n +2 | wc -l)
N_DEDUP=$(gzip -dc "$OUT_DEDUP" | tail -n +2 | wc -l)
N_PASS=$(gzip -dc "$OUT_PASS" | tail -n +2 | wc -l)
echo ""
echo "  sorted: $N_SORTED  dedup: $N_DEDUP  PASS: $N_PASS"
echo "=== done ==="
