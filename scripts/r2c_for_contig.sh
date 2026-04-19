#!/usr/bin/env bash
# r2c_for_contig.sh
#
# Given a svaba contig name, reconstruct the per-contig read-to-contig
# (r2c) alignment BAM for IGV. Uses the `bz:Z` aux tag that svaba writes
# on every r2c-aligned read (cname-keyed) to pick the right reads out
# of corrected.bam.
#
# Usage:
#   r2c_for_contig.sh CNAME CONTIGS_BAM CORRECTED_BAM [OUT_DIR]
#
# Example:
#   r2c_for_contig.sh c_chr12_10180001_10181001_0C \
#     T2_fermi.contigs.bam T2_fermi.corrected.bam \
#     r2c_debug/
#
# Output (in OUT_DIR):
#   contig.fa     single-contig reference FASTA (+ BWA & faidx indices)
#   reads.fq      corrected reads tagged bz:Z:<CNAME>
#   r2c.bam(.bai) sorted+indexed read-to-contig BAM — load this in IGV
#                 alongside contig.fa
#
# Tag semantics (v3):
#   bz:Z = cname list — every contig this read r2c-aligned to. This is
#          the right key for a per-contig read pull.
#   bi:Z = bp_id list — the specific variant rows this read supports
#          as ALT (matches r2c.txt.gz split_bps/disc_bps and bps.txt.gz
#          col 52). Pre-v3 this carried cnames; it no longer does.
#          If you want the ALT-supporter subset, set TAG=bi and pass
#          a bp_id as the first argument instead of a cname.

set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "usage: $0 CNAME CONTIGS_BAM CORRECTED_BAM [OUT_DIR]" >&2
  exit 2
fi

CNAME=$1
CONTIGS_BAM=$2
CORRECTED_BAM=$3
OUT_DIR=${4:-r2c_${CNAME}}

for t in samtools bwa awk; do
  command -v "$t" >/dev/null || { echo "ERROR: $t not on PATH" >&2; exit 1; }
done

mkdir -p "$OUT_DIR"
cd "$OUT_DIR"

# ---- 1. Pull the contig sequence out of contigs.bam -----------------------
# svaba writes each contig as a single BAM record whose QNAME == cname. Take
# SEQ as the reference. If the contig is reverse-complemented in the BAM
# (FLAG & 16), revcomp it back so the coordinate system matches cname.
# NB: use arithmetic bit extraction (int(flag/16)%2) instead of gawk's
# `and(flag,16)` so this works with BSD awk on macOS as well as gawk.
samtools view "../$CONTIGS_BAM" 2>/dev/null \
  | awk -v c="$CNAME" '
      BEGIN { found=0 }
      $1==c && !found {
        flag=$2+0; seq=$10
        if (int(flag/16) % 2 == 1) {
          # revcomp
          n=length(seq); out=""
          for (i=n;i>=1;i--) {
            b=substr(seq,i,1)
            if      (b=="A") out=out"T"
            else if (b=="T") out=out"A"
            else if (b=="C") out=out"G"
            else if (b=="G") out=out"C"
            else if (b=="a") out=out"t"
            else if (b=="t") out=out"a"
            else if (b=="c") out=out"g"
            else if (b=="g") out=out"c"
            else             out=out"N"
          }
          seq=out
        }
        print ">"c; print seq
        found=1
      }
      END { if (!found) exit 1 }
    ' > contig.fa || {
      echo "ERROR: contig $CNAME not found in $CONTIGS_BAM" >&2
      exit 1
    }

echo "contig length: $(awk 'NR==2{print length($0)}' contig.fa) bp"

# ---- 2. Build mini BWA index + faidx --------------------------------------
bwa index contig.fa 2> bwa_index.log
samtools faidx contig.fa

# ---- 3. Pull corrected reads tagged for this contig (boundary-aware) ------
# Default TAG=bz selects by cname (all r2c'd reads). TAG=bi selects by bp_id
# (ALT-supporter subset only) — in that mode, the CNAME positional arg is
# treated as a bp_id instead. See the v3 tag-semantics note in the header.
TAG=${TAG:-bz}
samtools view -h "../$CORRECTED_BAM" \
  | awk -v c="$CNAME" -v tag="$TAG" 'BEGIN{OFS="\t"; pfx="^"tag":Z:"}
      /^@/ { print; next }
      {
        for (i=12;i<=NF;i++) if ($i ~ pfx) {
          v=substr($i,6); n=split(v,a,",")
          for (j=1;j<=n;j++) if (a[j]==c) { print; next }
        }
      }' \
  | samtools fastq -@2 -n - > reads.fq 2> /dev/null

n_reads=$(( $(wc -l < reads.fq) / 4 ))
echo "reads tagged ${TAG}:Z:${CNAME}: $n_reads"

if [[ "$n_reads" -eq 0 ]]; then
  echo "ERROR: no reads tagged for $CNAME in $CORRECTED_BAM" >&2
  exit 1
fi

# ---- 4. Re-align reads to the single-contig reference ---------------------
# -k 15 lowers the minimum seed to pick up more partial matches on a small
# reference; -a outputs secondary hits so you can spot reads that map to
# multiple places on the contig (e.g. tandem repeats). Drop both if you
# want svaba-equivalent alignments.
bwa mem -t 4 -k 15 -a contig.fa reads.fq 2> bwa_mem.log \
  | samtools sort -@4 -o r2c.bam -
samtools index r2c.bam

echo
echo "Done. IGV: load contig.fa as Genome, then r2c.bam as track."
echo "Files in $PWD:"
ls -la r2c.bam r2c.bam.bai contig.fa contig.fa.fai reads.fq
