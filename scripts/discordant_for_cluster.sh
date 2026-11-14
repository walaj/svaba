#!/usr/bin/env bash
# discordant_for_cluster.sh
#
# Pull every read pair belonging to one svaba discordant cluster into a small
# sorted+indexed BAM for IGV, and (optionally) print the cluster's row from
# discordant.txt.gz so you can see the tumor/normal read split that drove the
# call. Handy for "why is this discordant cluster somatic when the normal has
# discordant reads?" debugging.
#
# The cluster id is the SAME string in all three places (svaba >= the
# unified-id change):
#   * bps.txt.gz  contig column  (for a DSCRD row, it replaces the contig)
#   * discordant.txt.gz  `id` column
#   * discordant.bam  `DC:Z:` read aux tag   (only present with --dump-reads)
# e.g.  FR_chr1_12345___chr5_67890
#
# Reads keep their source-prefixed QNAME (t001.../n001...), so tumor vs normal
# is visible in IGV by read name. The reads are at real genome coordinates —
# load them against your normal genome reference (hg38), NOT a contig.
#
# This is a samtools-only equivalent of:
#   svaba extract-pairs -i DISCORDANT_BAM -o OUT.bam --cluster CLUSTER_ID
# Use the subcommand if you have a current svaba build; use this script if you
# only have samtools, or you also want the discordant.txt.gz context line.
#
# Usage:
#   discordant_for_cluster.sh CLUSTER_ID DISCORDANT_BAM [DISCORDANT_TXT_GZ] [OUT_DIR]
#
# Example:
#   discordant_for_cluster.sh 'FR_chr1_12345___chr5_67890' \
#     sample.discordant.bam sample.discordant.txt.gz disc_debug/

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "usage: $0 CLUSTER_ID DISCORDANT_BAM [DISCORDANT_TXT_GZ] [OUT_DIR]" >&2
  exit 2
fi

CLUSTER_ID=$1
DISC_BAM=$2
DISC_TXT=${3:-}
SAN=$(printf '%s' "$CLUSTER_ID" | tr -c 'A-Za-z0-9._-' '_')
OUT_DIR=${4:-disc_${SAN}}

command -v samtools >/dev/null || { echo "ERROR: samtools not on PATH" >&2; exit 1; }
[[ -f "$DISC_BAM" ]] || { echo "ERROR: discordant BAM not found: $DISC_BAM" >&2; exit 1; }

mkdir -p "$OUT_DIR"

# ---- 0. (optional) show the cluster's discordant.txt.gz row ----------------
# Columns: chr1 pos1 strand1 chr2 pos2 strand2 tcount ncount mapq1 mapq2
#          nm1 nm2 cname id valid. tcount/ncount = tumor/normal read counts.
if [[ -n "$DISC_TXT" && -f "$DISC_TXT" ]]; then
  echo "=== discordant.txt.gz rows for $CLUSTER_ID ==="
  # header (line 1) + any row whose id column matches exactly
  { gzip -dc "$DISC_TXT" 2>/dev/null || zcat "$DISC_TXT"; } \
    | awk -F'\t' -v id="$CLUSTER_ID" 'NR==1 || $14==id'
  echo
fi

# ---- 1. pass 1: collect QNAMEs whose DC:Z tag contains this cluster id ------
# Literal substring match via index() (the id has no regex metachars, and a
# cluster id can contain 'd' from a "*_decoy" contig so we don't tokenize the
# "d"-joined DC list). The "___" inside every id makes prefix collisions
# (e.g. _100 vs _1000) practically impossible.
QN="$OUT_DIR/qnames.txt"
samtools view "$DISC_BAM" \
  | awk -v id="$CLUSTER_ID" '{
      for (i=12;i<=NF;i++)
        if (substr($i,1,5)=="DC:Z:" && index($i, id)) { print $1; break }
    }' \
  | sort -u > "$QN"

n_qn=$(wc -l < "$QN" | awk '{print $1}')
echo "qnames carrying DC:Z containing $CLUSTER_ID: $n_qn"
if [[ "$n_qn" -eq 0 ]]; then
  echo "ERROR: no reads in $DISC_BAM carry that cluster id." >&2
  echo "       (was svaba run with --dump-reads? is the id spelled exactly as" >&2
  echo "        in the bps.txt contig column / discordant.txt.gz id column?)" >&2
  exit 1
fi

# ---- 2. pass 2: extract all alignments for those qnames, sort, index --------
# samtools view -N <file> keeps records whose read name is listed (both mates
# + any supplementary share the QNAME, so the pairs stay together).
OUT_BAM="$OUT_DIR/cluster.bam"
samtools view -h -N "$QN" "$DISC_BAM" \
  | samtools sort -@4 -o "$OUT_BAM" -
samtools index "$OUT_BAM"

n_rec=$(samtools view -c "$OUT_BAM")
echo
echo "Done. $n_rec records -> $OUT_BAM"
echo "IGV: load your genome reference (hg38), then $OUT_BAM as a track."
echo "     Tumor vs normal reads are distinguished by the QNAME prefix"
echo "     (t001.../n001...); colour/group by read name in IGV."
ls -la "$OUT_BAM" "$OUT_BAM.bai"
