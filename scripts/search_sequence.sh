#!/bin/bash
# ============================================================
# Search a BAM for reads matching any of the given sequences
# (or their reverse complements) and write matches to a new BAM.
#
# Usage:
#   ./bam_grep.sh <input.bam> <output_prefix> <seq1> [seq2 ...]
#
# Example:
#   ./bam_grep.sh tumor.bam egfr_vIII ACTGACTGACTG TGCATGCATGCA
#   -> produces egfr_vIII.bam + egfr_vIII.bam.bai
# ============================================================

set -euo pipefail

# ---- Args ----
if [ $# -lt 3 ]; then
    echo "Usage: $0 <input.bam> <output_prefix> <seq1> [seq2 ...]" >&2
    echo "  Searches for each sequence and its reverse complement." >&2
    echo "  Output: <output_prefix>.bam + .bam.bai" >&2
    exit 1
fi

INPUT_BAM="$1"
OUT_PREFIX="$2"
shift 2
SEQS=("$@")

OUT_BAM="${OUT_PREFIX}.bam"
THREADS=${THREADS:-4}

# ---- Sanity checks ----
if [ ! -f "$INPUT_BAM" ]; then
    echo "ERROR: Input BAM not found: $INPUT_BAM" >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found in PATH" >&2
    exit 1
fi

# ---- Build the combined regex: seq1|revcomp1|seq2|revcomp2|... ----
revcomp() {
    # Uppercase, complement, reverse
    echo "$1" | tr 'acgtnACGTN' 'tgcanTGCAN' | rev
}

PATTERNS=()
echo "Searching for sequences (and reverse complements):" >&2
for seq in "${SEQS[@]}"; do
    # Uppercase and strip whitespace
    seq_up=$(echo "$seq" | tr 'a-z' 'A-Z' | tr -d '[:space:]')
    rc=$(revcomp "$seq_up")
    PATTERNS+=("$seq_up" "$rc")
    echo "  $seq_up  (rc: $rc)" >&2
done

# Join patterns with | for egrep
PATTERN=$(IFS='|'; echo "${PATTERNS[*]}")

echo "" >&2
echo "Input:   $INPUT_BAM" >&2
echo "Output:  $OUT_BAM" >&2
echo "Threads: $THREADS" >&2
echo "" >&2

# ---- Filter ----
# samtools view streams SAM (header + reads), awk keeps headers and any read
# whose sequence column (field 10) matches the pattern, samtools view -b
# re-encodes BAM.
samtools view -@ "$THREADS" -h "$INPUT_BAM" \
    | awk -v pat="$PATTERN" '
        /^@/ { print; next }
        $10 ~ pat { print }
      ' \
    | samtools view -@ "$THREADS" -b -o "$OUT_BAM" -

# ---- Sort (safer - preserves coordinate order for indexing) ----
# If input was already sorted, output is too (we only filter), but re-sort
# defensively in case the original had unsorted regions.
samtools sort -@ "$THREADS" -o "${OUT_BAM}.sorted" "$OUT_BAM"
mv "${OUT_BAM}.sorted" "$OUT_BAM"

# ---- Index ----
samtools index -@ "$THREADS" "$OUT_BAM"

# ---- Report ----
N_READS=$(samtools view -c "$OUT_BAM")
echo "Done. Matched reads: $N_READS" >&2
echo "Output: $OUT_BAM + ${OUT_BAM}.bai" >&2
