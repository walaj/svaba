#!/usr/bin/env bash
#
# mix_tumor_normal.sh -- in-silico tumor impurity (admixture).
#
# Builds an IMPURE tumor BAM by mixing simulated-tumor reads with normal reads
# at a target purity and depth:
#     impure_tumor = purity * (tumor reads) + (1-purity) * (normal reads)
# A somatic variant that is heterozygous in the pure tumor (VAF 0.5) appears at
# VAF ~ 0.5*purity in the mix. The unmixed normal BAM remains the matched normal
# for the svaba T/N pair.
#
# Usage:
#   mix_tumor_normal.sh -t TUMOR.bam -n NORMAL.bam -p PURITY -d TARGET_DEPTH -o OUT.bam \
#                       [--tumor-cov X] [--normal-cov Y] [--readlen 150] [--seed 42] [--threads 8]
#
# If --tumor-cov / --normal-cov are omitted they are estimated from the BAM via
# `samtools idxstats` (mapped_reads * readlen / genome_len).
set -euo pipefail
SAMTOOLS="${SAMTOOLS:-samtools}"
READLEN=150; SEED=42; THREADS=8; TUMOR_COV=""; NORMAL_COV=""
T=""; N=""; P=""; D=""; OUT=""

while [ $# -gt 0 ]; do
  case "$1" in
    -t) T="$2"; shift 2;;
    -n) N="$2"; shift 2;;
    -p) P="$2"; shift 2;;
    -d) D="$2"; shift 2;;
    -o) OUT="$2"; shift 2;;
    --tumor-cov) TUMOR_COV="$2"; shift 2;;
    --normal-cov) NORMAL_COV="$2"; shift 2;;
    --readlen) READLEN="$2"; shift 2;;
    --seed) SEED="$2"; shift 2;;
    --threads) THREADS="$2"; shift 2;;
    -h|--help) grep '^#' "$0" | sed 's/^# \{0,1\}//'; exit 0;;
    *) echo "unknown arg: $1"; exit 1;;
  esac
done
[ -n "$T" ] && [ -n "$N" ] && [ -n "$P" ] && [ -n "$D" ] && [ -n "$OUT" ] || {
  echo "ERROR: need -t -n -p -d -o (see --help)"; exit 1; }
log(){ printf '[mix %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

est_cov(){  # $1=bam -> mean depth over mapped contigs
  "$SAMTOOLS" idxstats "$1" | awk -v rl="$READLEN" '
    $1!="*"{len+=$2; mapped+=$3} END{ if(len>0) printf "%.3f", mapped*rl/len; else print 0 }'
}
[ -z "$TUMOR_COV" ]  && { log "estimating tumor coverage ..."; TUMOR_COV=$(est_cov "$T"); }
[ -z "$NORMAL_COV" ] && { log "estimating normal coverage ..."; NORMAL_COV=$(est_cov "$N"); }
log "tumor_cov=$TUMOR_COV  normal_cov=$NORMAL_COV  purity=$P  target_depth=$D"

# fraction of each source BAM to keep
FT=$(awk -v p="$P" -v d="$D" -v c="$TUMOR_COV"  'BEGIN{printf "%.5f", (p*d)/c}')
FN=$(awk -v p="$P" -v d="$D" -v c="$NORMAL_COV" 'BEGIN{printf "%.5f", ((1-p)*d)/c}')
log "keep tumor frac=$FT  normal frac=$FN"
awk -v f="$FT" 'BEGIN{exit !(f<=1.0)}' || log "WARN: need more tumor depth (frac>1, clamping to 1.0)"
awk -v f="$FN" 'BEGIN{exit !(f<=1.0)}' || log "WARN: need more normal depth (frac>1, clamping to 1.0)"
FT=$(awk -v f="$FT" 'BEGIN{print (f>1)?1.0:f}')
FN=$(awk -v f="$FN" 'BEGIN{print (f>1)?1.0:f}')

TMP=$(mktemp -d "${OUT%.bam}.mixtmp.XXXX")
trap 'rm -rf "$TMP"' EXIT
RG='@RG\tID:impure_tumor\tSM:tumor\tLB:impure\tPL:ILLUMINA'

sub(){  # $1=bam $2=frac $3=out  -- subsample + force a single shared RG
  local f="$2"
  if awk -v f="$f" 'BEGIN{exit !(f>=0.9999)}'; then
    "$SAMTOOLS" addreplacerg -@ "$THREADS" -r "$RG" -o "$3" "$1"
  else
    "$SAMTOOLS" view -@ "$THREADS" -b -s "${SEED}${f#0}" "$1" \
      | "$SAMTOOLS" addreplacerg -@ "$THREADS" -r "$RG" -o "$3" -
  fi
}
log "subsampling tumor ..."  ; sub "$T" "$FT" "$TMP/t.bam"
log "subsampling normal ..." ; sub "$N" "$FN" "$TMP/n.bam"
log "merging -> $OUT ..."
"$SAMTOOLS" merge -@ "$THREADS" -c -f -o - "$TMP/t.bam" "$TMP/n.bam" \
  | "$SAMTOOLS" sort -@ "$THREADS" -o "$OUT" -
"$SAMTOOLS" index -@ "$THREADS" "$OUT"
log "DONE. impure tumor @ purity=$P depth~=$D : $OUT"
log "  pair with the UNMIXED normal: svaba run -t $OUT -n $N -G \$REF ..."
