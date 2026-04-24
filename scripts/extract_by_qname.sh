#!/bin/bash
# ============================================================
# Extract all reads from a target BAM whose QNAMEs appear in
# a query BAM.
#
# Both inputs can be coordinate-sorted. We extract QNAMEs from
# the query BAM into a hash set, then stream the target BAM and
# keep any record whose QNAME is in the set. No requirement that
# the BAMs are sorted the same way, or at all.
#
# Usage:
#   ./extract_by_qname.sh <query.bam> <target.bam> <output.bam>
#
# Example:
#   ./extract_by_qname.sh matches.bam /mnt/ssd/tumor.recal.bam out.bam
# ============================================================

set -euo pipefail

if [ $# -ne 3 ]; then
    echo "Usage: $0 <query.bam> <target.bam> <output.bam>" >&2
    echo "  Extracts from target.bam all reads whose QNAME appears" >&2
    echo "  in query.bam. Writes sorted, indexed BAM to output.bam." >&2
    exit 1
fi

QUERY_BAM="$1"
TARGET_BAM="$2"
OUT_BAM="$3"
THREADS=${THREADS:-8}

# ---- Sanity checks ----
for f in "$QUERY_BAM" "$TARGET_BAM"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: BAM not found: $f" >&2
        exit 1
    fi
done

if ! command -v samtools >/dev/null 2>&1; then
    echo "ERROR: samtools not found in PATH" >&2
    exit 1
fi

# ---- Use a scratch dir for intermediate files ----
TMPDIR=$(mktemp -d -t extract_qname.XXXXXX)
trap "rm -rf $TMPDIR" EXIT
QNAME_LIST="$TMPDIR/qnames.txt"

echo "Temp dir: $TMPDIR" >&2

# ---- Step 1: extract unique QNAMEs from query BAM ----
echo "[1/3] Extracting QNAMEs from $QUERY_BAM..." >&2
samtools view -@ "$THREADS" "$QUERY_BAM" \
    | awk '{print $1}' \
    | sort -u \
    > "$QNAME_LIST"

N_QNAMES=$(wc -l < "$QNAME_LIST")
echo "  Found $N_QNAMES unique QNAMEs" >&2

if [ "$N_QNAMES" -eq 0 ]; then
    echo "ERROR: No QNAMEs found in query BAM" >&2
    exit 1
fi

# ---- Step 2: stream target BAM, filter by QNAME ----
echo "[2/3] Filtering $TARGET_BAM..." >&2
# Use awk with a hash lookup - O(1) per read, linear in target BAM size
samtools view -@ "$THREADS" -h "$TARGET_BAM" \
    | awk -v list="$QNAME_LIST" '
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

# ---- Step 3: sort and index ----
echo "[3/3] Sorting and indexing output..." >&2
samtools sort -@ "$THREADS" -o "$OUT_BAM" "$TMPDIR/unsorted.bam"
samtools index -@ "$THREADS" "$OUT_BAM"

# ---- Report ----
N_OUT=$(samtools view -c "$OUT_BAM")
echo "" >&2
echo "Done." >&2
echo "  Query QNAMEs:       $N_QNAMES" >&2
echo "  Output reads:       $N_OUT" >&2
echo "  (Output includes all alignments per QNAME - primary, secondary," >&2
echo "   supplementary, and both mates of pairs.)" >&2
echo "  Output: $OUT_BAM + ${OUT_BAM%.bam}.bam.bai" >&2
