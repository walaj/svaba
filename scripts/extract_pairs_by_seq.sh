#!/bin/bash
# ============================================================
# Extract read pairs from a BAM where either mate contains any
# of the given sequences (or their reverse complements).
#
# Two-pass approach on the target BAM:
#   Pass 1: stream through BAM, collect QNAMEs of reads whose
#           SEQ matches any query sequence or its reverse complement
#   Pass 2: stream through BAM again, extract ALL reads (both mates,
#           plus supplementary/secondary) whose QNAME was collected
#
# Usage:
#   ./extract_pairs_by_seq.sh <target.bam> <output.bam> <seq1> [seq2 ...]
#
# Example:
#   ./extract_pairs_by_seq.sh tumor.bam fusion_reads.bam \
#       ACGTACGTACGT TTACGATCGATC
# ============================================================

set -euo pipefail

# ---- Args ----
if [ $# -lt 3 ]; then
    echo "Usage: $0 <target.bam> <output.bam> <seq1> [seq2 ...]" >&2
    echo "  Extracts read pairs where either mate matches any sequence" >&2
    echo "  (or its reverse complement)." >&2
    echo "  Output: sorted, indexed BAM + BAI." >&2
    exit 1
fi

TARGET_BAM="$1"
OUT_BAM="$2"
shift 2
SEQS=("$@")

THREADS=${THREADS:-8}
AWK_BIN=$(command -v mawk || command -v awk)

# ---- Sanity checks ----
if [ ! -f "$TARGET_BAM" ]; then
    echo "ERROR: Target BAM not found: $TARGET_BAM" >&2
    exit 1
fi

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found in PATH" >&2
    exit 1
fi

# ---- Scratch dir ----
TMPDIR=$(mktemp -d -t extract_pairs.XXXXXX)
trap "rm -rf $TMPDIR" EXIT
QNAME_LIST="$TMPDIR/qnames.txt"

echo "Target: $TARGET_BAM" >&2
echo "Output: $OUT_BAM" >&2
echo "Temp:   $TMPDIR" >&2
echo "Awk:    $AWK_BIN" >&2
echo "Threads: $THREADS" >&2
echo "" >&2

# ---- Build combined regex: seq1|revcomp1|seq2|revcomp2|... ----
revcomp() {
    echo "$1" | tr 'acgtnACGTN' 'tgcanTGCAN' | rev
}

PATTERNS=()
echo "Query sequences (with reverse complements):" >&2
for seq in "${SEQS[@]}"; do
    seq_up=$(echo "$seq" | tr 'a-z' 'A-Z' | tr -d '[:space:]')
    rc=$(revcomp "$seq_up")
    PATTERNS+=("$seq_up" "$rc")
    echo "  $seq_up  (rc: $rc)" >&2
done
PATTERN=$(IFS='|'; echo "${PATTERNS[*]}")
echo "" >&2

# ---- PASS 1: find QNAMEs whose SEQ matches ----
echo "[Pass 1/2] Scanning $TARGET_BAM for matching sequences..." >&2
samtools view -@ "$THREADS" "$TARGET_BAM" \
    | "$AWK_BIN" -v pat="$PATTERN" '$10 ~ pat {print $1}' \
    | sort -u \
    > "$QNAME_LIST"

N_QNAMES=$(wc -l < "$QNAME_LIST")
echo "  Matching QNAMEs: $N_QNAMES" >&2

if [ "$N_QNAMES" -eq 0 ]; then
    echo "No matches found. Exiting without creating output." >&2
    exit 0
fi

# ---- PASS 2: extract all alignments for those QNAMEs (both mates) ----
echo "[Pass 2/2] Extracting all alignments for matching QNAMEs..." >&2
samtools view -@ "$THREADS" -h "$TARGET_BAM" \
    | "$AWK_BIN" -v list="$QNAME_LIST" '
        BEGIN {
            while ((getline line < list) > 0) {
                qnames[line] = 1
            }
            close(list)
        }
        /^@/ { print; next }
        ($1 in qnames) { print }
      ' \
    | samtools view -@ "$THREADS" -b -o "$TMPDIR/unsorted.bam" -

# ---- Sort and index ----
echo "Sorting and indexing..." >&2
samtools sort -@ "$THREADS" -o "$OUT_BAM" "$TMPDIR/unsorted.bam"
samtools index -@ "$THREADS" "$OUT_BAM"

# ---- Report ----
N_OUT=$(samtools view -c "$OUT_BAM")
N_PRIMARY=$(samtools view -c -F 2304 "$OUT_BAM")  # exclude secondary + supplementary

echo "" >&2
echo "Done." >&2
echo "  Matching QNAMEs (read pairs):  $N_QNAMES" >&2
echo "  Total alignments in output:    $N_OUT" >&2
echo "  Primary alignments (both mates): $N_PRIMARY" >&2
echo "  Output: $OUT_BAM + ${OUT_BAM%.bam}.bam.bai" >&2
