#!/usr/bin/env bash
#
# run_sim.sh -- end-to-end svaba SV simulation:
#   1. generate non-overlapping SVs + diploid donor genome + truth set
#   2. ART -> paired FASTQ from each haplotype (half coverage each)
#   3. BWA-MEM align back to the reference -> sorted/indexed tumor BAM @ MAXCOV
#   4. downsample to each requested coverage (samtools view -s)
#   5. (optional) simulate a matched clean normal the same way
#
# Outputs (in $OUTDIR):
#   tumor.hap1.fa tumor.hap2.fa          donor haplotypes
#   truth.vcf truth.bedpe truth.sv_table.tsv
#   tumor.<cov>x.bam(.bai)               aligned tumor at each coverage
#   normal.<cov>x.bam(.bai)              (only if SIM_NORMAL=1)
#
# Tumor/normal pair for svaba:
#   svaba run -t tumor.<cov>x.bam -n <normal.bam> -G $REF ...
# Impurity: see mix_tumor_normal.sh.
#
# Config: edit sim_config.sh or export overrides before calling.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/sim_config.sh"

log(){ printf '[run_sim %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }
mkdir -p "$OUTDIR"

# sanity
command -v "$ART" >/dev/null || { echo "ERROR: art_illumina not found ($ART)"; exit 1; }
command -v "$BWA" >/dev/null || { echo "ERROR: bwa not found"; exit 1; }
command -v "$SAMTOOLS" >/dev/null || { echo "ERROR: samtools not found"; exit 1; }
[ -f "$REF.fai" ]  || { echo "ERROR: $REF.fai missing (samtools faidx $REF)"; exit 1; }
[ -f "$REF.bwt" ]  || { echo "ERROR: $REF.bwt missing (bwa index $REF)"; exit 1; }

HALFCOV=$(awk -v c="$MAXCOV" 'BEGIN{printf "%.4f", c/2.0}')

# ---------------------------------------------------------------------------
# 1. donor genome + truth
# ---------------------------------------------------------------------------
if [ -f "$OUTDIR/truth.bedpe" ] && [ -f "$OUTDIR/tumor.hap1.fa" ]; then
  log "donor genome + truth already present -> skipping generation"
else
  log "generating SVs + donor genome (seed=$SEED, chroms=$CHROMS) ..."
  avoid_arg=()
  [ -n "$AVOID_BED" ] && avoid_arg=(--avoid-bed "$AVOID_BED")
  # bash 3.2 (macOS): guard empty-array expansion under `set -u`
  "$PYTHON" "$SIM_PY" --ref "$REF" --out-dir "$OUTDIR" --seed "$SEED" \
      --chroms "$CHROMS" ${avoid_arg[@]+"${avoid_arg[@]}"} $GEN_ARGS
fi

# ---------------------------------------------------------------------------
# helper: simulate one diploid sample (two hap FASTAs) at TOTAL coverage and
# align to REF. $1=sample tag (tumor/normal) $2=hap1.fa $3=hap2.fa $4=total cov
# $5=output bam
# ---------------------------------------------------------------------------
simulate_and_align(){
  local tag="$1" h1="$2" h2="$3" cov="$4" outbam="$5"
  local half; half=$(awk -v c="$cov" 'BEGIN{printf "%.4f", c/2.0}')
  local wd="$OUTDIR/fastq_$tag"; mkdir -p "$wd"

  for h in 1 2; do
    local fa; [ "$h" = 1 ] && fa="$h1" || fa="$h2"
    if [ -f "$wd/${tag}_h${h}_1.fq" ] || [ -f "$wd/${tag}_h${h}_1.fq.gz" ]; then
      log "$tag hap$h FASTQ exists -> skip ART"
    else
      log "ART: $tag hap$h @ ${half}x ($(basename "$fa")) ..."
      "$ART" -ss "$ART_PROFILE" -na -i "$fa" -p -l "$READLEN" -f "$half" \
             -m "$FRAG_MEAN" -s "$FRAG_SD" -rs "$((SEED + h))" \
             -o "$wd/${tag}_h${h}_" >/dev/null
    fi
  done

  log "BWA-MEM align $tag @ ${cov}x -> $(basename "$outbam") ..."
  local rg="@RG\tID:${tag}\tSM:${tag}\tLB:${tag}\tPL:ILLUMINA"
  # concat both haplotypes' reads on the fly
  "$BWA" mem -t "$THREADS" -R "$rg" "$REF" \
      <(cat "$wd/${tag}_h1_1.fq" "$wd/${tag}_h2_1.fq") \
      <(cat "$wd/${tag}_h1_2.fq" "$wd/${tag}_h2_2.fq") 2>"$OUTDIR/bwa_$tag.log" \
    | "$SAMTOOLS" sort -@ "$THREADS" -o "$outbam" -
  "$SAMTOOLS" index -@ "$THREADS" "$outbam"
}

# ---------------------------------------------------------------------------
# 2-4. tumor: simulate @ MAXCOV, then downsample
# ---------------------------------------------------------------------------
TUMOR_MAX_BAM="$OUTDIR/tumor.${MAXCOV}x.bam"
if [ -f "$TUMOR_MAX_BAM" ]; then
  log "tumor @ ${MAXCOV}x exists -> skip"
else
  simulate_and_align tumor "$OUTDIR/tumor.hap1.fa" "$OUTDIR/tumor.hap2.fa" "$MAXCOV" "$TUMOR_MAX_BAM"
fi

for cov in $COVERAGES; do
  out="$OUTDIR/tumor.${cov}x.bam"
  [ "$cov" = "$MAXCOV" ] && continue
  if [ -f "$out" ]; then log "tumor @ ${cov}x exists -> skip"; continue; fi
  frac=$(awk -v t="$cov" -v m="$MAXCOV" 'BEGIN{printf "%.4f", t/m}')
  awk -v f="$frac" 'BEGIN{exit !(f<=1.0)}' || { log "WARN: ${cov}x > MAXCOV; skipping"; continue; }
  log "downsample tumor ${MAXCOV}x -> ${cov}x (keep $frac) ..."
  "$SAMTOOLS" view -@ "$THREADS" -b -s "${SEED}.${frac#0.}" "$TUMOR_MAX_BAM" -o "$out"
  "$SAMTOOLS" index -@ "$THREADS" "$out"
done

# ---------------------------------------------------------------------------
# 5. normal (optional): simulate a clean normal from the reference itself.
#    Default is to use a real NORMAL_BAM instead (no simulation).
# ---------------------------------------------------------------------------
if [ "$SIM_NORMAL" = "1" ]; then
  # If germline donor exists use it; otherwise simulate from the bare reference
  # (split into two identical haplotypes => clean diploid normal).
  if [ -f "$OUTDIR/normal.hap1.fa" ]; then
    NH1="$OUTDIR/normal.hap1.fa"; NH2="$OUTDIR/normal.hap2.fa"
  else
    log "no germline donor; building clean normal haplotypes from REF subset ..."
    # extract chosen chroms from REF once for a clean diploid normal
    if [ ! -f "$OUTDIR/ref_subset.fa" ]; then
      "$PYTHON" - "$REF" "$OUTDIR/ref_subset.fa" "$CHROMS" <<'PY'
import sys,subprocess
ref,out,chroms=sys.argv[1],sys.argv[2],sys.argv[3]
cs=[]
for tok in chroms.split(','):
    tok=tok.strip()
    if '-' in tok and tok.lower().startswith('chr') and tok[3:].split('-')[0].isdigit():
        lo,hi=tok[3:].split('-'); cs+=["chr%d"%i for i in range(int(lo),int(hi)+1)]
    elif tok: cs.append(tok)
with open(out,'w') as o:
    for c in cs:
        o.write(subprocess.run(["samtools","faidx",ref,c],capture_output=True,text=True,check=True).stdout)
PY
    fi
    NH1="$OUTDIR/ref_subset.fa"; NH2="$OUTDIR/ref_subset.fa"
  fi
  NORMAL_MAX_BAM="$OUTDIR/normal.${NORMAL_COV}x.bam"
  if [ -f "$NORMAL_MAX_BAM" ]; then log "normal @ ${NORMAL_COV}x exists -> skip"; else
    simulate_and_align normal "$NH1" "$NH2" "$NORMAL_COV" "$NORMAL_MAX_BAM"
  fi
  log "simulated normal: $NORMAL_MAX_BAM"
else
  log "using real normal: $NORMAL_BAM (set SIM_NORMAL=1 to simulate one instead)"
fi

log "DONE. Tumor BAMs:"; ls -1 "$OUTDIR"/tumor.*x.bam 2>/dev/null | sed 's/^/   /' >&2
log "Truth: $OUTDIR/truth.bedpe  |  Example svaba T/N run:"
n="${NORMAL_BAM}"; [ "$SIM_NORMAL" = "1" ] && n="$OUTDIR/normal.${NORMAL_COV}x.bam"
log "   svaba run -t $OUTDIR/tumor.${MAXCOV}x.bam -n $n -G $REF -a sim -p $THREADS"
