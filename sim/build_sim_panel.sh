#!/usr/bin/env bash
#
# build_sim_panel.sh -- build a reproducible, self-contained simulated tumor panel:
#   * many SVs on chr22 (default) -> diploid donor genome + ground truth
#   * ART -> BWA-MEM simulated tumor reads
#   * in-silico tumor IMPURITY by admixing a REAL whole-genome normal
#     (downsampled), at several coverages x purities
#   * one new timestamped+seeded output folder with the indexed tumor BAMs,
#     all ground-truth files (VCF/BEDPE/TSV), and a README documenting exactly
#     how it was made (incl. the seed) so it can be regenerated.
#
# Reuses simulate_sv_genome.py + run_sim.sh + mix_tumor_normal.sh (this dir).
#
# Why a SECOND normal for contamination: using a different real normal from the
# SAME person as the matched normal means the tumor's normal-cell contamination
# shares the individual's germline WITHOUT reusing the matched normal's exact
# reads -- avoids circular/identical-read artifacts at germline sites.
#
# Run it (when ready):   ./build_sim_panel.sh
# Preview the plan only:  DRYRUN=1 ./build_sim_panel.sh
# Quick validation:       SMOKE=1 ./build_sim_panel.sh   (tiny counts/cov, fast)
#
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ----------------------------------------------------------------------------
# config (override any via environment)
# ----------------------------------------------------------------------------
: "${OUTROOT:=/Volumes/wala24T/sim}"          # each run -> a new subfolder here
: "${CACHE:=$OUTROOT/cache}"                   # reusable normal-pool cache
: "${REF:=$HOME/ref_genome/Homo_sapiens_assembly38.fasta}"
: "${THREADS:=10}"
: "${SEED:=42}"                               # recorded in the README
: "${TAG:=wgs}"

: "${CHROMS:=chr1-22,chrX,chrY}"              # whole genome (autosomes + X + Y)
: "${COVERAGES:=5 10 15 20}"                  # total (genome-wide) depth of each tumor BAM
: "${PURITIES:=0.5}"                          # tumor fraction; (1-p) = normal admixture
: "${NORMAL_DOWNSAMPLE:=0.5}"                 # randomly keep this frac of the real normal as the contamination pool

# contamination normal (a 2nd WGS normal from the same individual)
: "${CONTAM_NORMAL_BAM:=/Volumes/wala24T/ss/genomics/genomics-bulk/2022.12.16/DNA/2022.12.16.dna.personalis.WGS/aligned_bams/recalibrated/NP_WGS_Normal_DNA/NP_WGS_Normal_DNA.recal.bam}"
# matched normal for the svaba T/N pair (the OTHER normal). README-only; optional.
: "${MATCHED_NORMAL_BAM:=$HOME/Downloads/svaba_compare/blood.recal.bam}"

# ART params (forwarded to run_sim.sh)
: "${READLEN:=150}"; : "${FRAG_MEAN:=400}"; : "${FRAG_SD:=50}"; : "${ART_PROFILE:=HS25}"

# how many SVs: thousands across the whole genome. ~2.9 Gb callable with
# min-gap 50 kb keeps events independent (>= svaba's assembly window). Includes
# inter-chromosomal translocations (n-tra) since >1 chrom is present. ~5,200 SVs;
# scale --n-* up freely. Sizes use the generator defaults (DEL 50bp-3Mb,
# DUP/INV 200bp-1Mb, TRA segment 1kb-500kb).
: "${GEN_ARGS:=--n-del-small 1200 --n-ins-small 1200 --n-del 1000 --n-dup 700 --n-inv 700 --n-tra 400 --min-gap 50000 --hom-frac 0.3 --germline-frac 0.0}"

: "${DRYRUN:=0}"
: "${SMOKE:=0}"
if [ "$SMOKE" = "1" ]; then           # fast end-to-end validation (2 small chroms)
  TAG="smoke"; CHROMS="chr21,chr22"; COVERAGES="5 10"; PURITIES="0.5"
  GEN_ARGS="--n-del-small 60 --n-ins-small 60 --n-del 40 --n-dup 30 --n-inv 30 --n-tra 10 --min-gap 15000 --del-max 1000000 --dup-max 500000 --inv-max 500000"
fi

: "${SAMTOOLS:=samtools}"; : "${PYTHON:=python3}"
log(){ printf '[panel %s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

# ----------------------------------------------------------------------------
# preflight
# ----------------------------------------------------------------------------
[ -f "$REF.fai" ] || { echo "ERROR: $REF.fai missing"; exit 1; }
[ -f "$REF.bwt" ] || { echo "ERROR: $REF.bwt missing (bwa index)"; exit 1; }
[ -f "$CONTAM_NORMAL_BAM" ] || { echo "ERROR: contamination normal not found: $CONTAM_NORMAL_BAM"; exit 1; }
[ -f "$CONTAM_NORMAL_BAM.bai" ] || log "WARN: no .bai for contamination normal (idxstats/subsample may be slow)"
mkdir -p "$OUTROOT" "$CACHE"

RUNDIR="$OUTROOT/${TAG}_$(date +%Y%m%d_%H%M%S)_seed${SEED}"
# tumor sim coverage: must cover p*C for the densest cell in the grid (+margin)
maxC=$(printf '%s\n' $COVERAGES | sort -n | tail -1)
maxP=$(printf '%s\n' $PURITIES  | sort -n | tail -1)
TUMOR_SIM_COV=$(awk -v c="$maxC" -v p="$maxP" 'BEGIN{printf "%d", (c*p)*1.15 + 2}')

# expand CHROMS (chr1-22 ranges) -> space list, for coverage-matching the
# contamination over the simulated chromosomes (genome-wide for the WGS default).
expand_chroms(){ local out=() tok i; for tok in $(printf '%s' "$1" | tr ',' ' '); do
  if [[ "$tok" =~ ^chr([0-9]+)-([0-9]+)$ ]]; then
    for i in $(seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"); do out+=("chr$i"); done
  else out+=("$tok"); fi; done; printf '%s ' "${out[@]}"; }
CHROM_LIST=$(expand_chroms "$CHROMS")

echo "============================================================" >&2
log "PLAN"
log "  output folder : $RUNDIR"
log "  reference     : $REF"
log "  SV chroms     : $CHROMS   (seed $SEED)"
log "  SV counts     : $GEN_ARGS"
log "  tumor sim cov : ${TUMOR_SIM_COV}x  (covers max p*C = ${maxP}*${maxC})"
log "  coverages     : $COVERAGES   purities: $PURITIES"
log "  contamination : $CONTAM_NORMAL_BAM (downsample ${NORMAL_DOWNSAMPLE})"
log "  threads       : $THREADS"
N=0; for C in $COVERAGES; do for P in $PURITIES; do N=$((N+1)); done; done
log "  -> $N impure tumor BAM(s): tumor.<cov>x.p<purity>.bam"
echo "============================================================" >&2
if [ "$DRYRUN" = "1" ]; then log "DRYRUN=1 -> stopping before any heavy work."; exit 0; fi
mkdir -p "$RUNDIR"

# ----------------------------------------------------------------------------
# 1+2. donor genome + truth + pure tumor reads (delegate to run_sim.sh)
# ----------------------------------------------------------------------------
log "generating SVs + donor genome + pure tumor @ ${TUMOR_SIM_COV}x ..."
OUTDIR="$RUNDIR" CHROMS="$CHROMS" SEED="$SEED" THREADS="$THREADS" REF="$REF" \
  MAXCOV="$TUMOR_SIM_COV" COVERAGES="$TUMOR_SIM_COV" SIM_NORMAL=0 \
  READLEN="$READLEN" FRAG_MEAN="$FRAG_MEAN" FRAG_SD="$FRAG_SD" ART_PROFILE="$ART_PROFILE" \
  GEN_ARGS="$GEN_ARGS" \
  "$HERE/run_sim.sh"
PURE_TUMOR="$RUNDIR/tumor.${TUMOR_SIM_COV}x.bam"
[ -f "$PURE_TUMOR" ] || { echo "ERROR: pure tumor not produced ($PURE_TUMOR)"; exit 1; }

# ----------------------------------------------------------------------------
# 3. contamination pool: randomly keep NORMAL_DOWNSAMPLE of the real normal
#    (cached + reused across runs; this is the literal "downsample 50% then mix")
# ----------------------------------------------------------------------------
POOL="$CACHE/$(basename "${CONTAM_NORMAL_BAM%.bam}").ds${NORMAL_DOWNSAMPLE}.seed${SEED}.bam"
if [ -f "$POOL" ] && [ -f "$POOL.bai" ]; then
  log "contamination pool exists -> reuse $(basename "$POOL")"
else
  log "downsampling contamination normal to ${NORMAL_DOWNSAMPLE} (one-time, large) ..."
  "$SAMTOOLS" view -@ "$THREADS" -b -s "${SEED}.${NORMAL_DOWNSAMPLE#0.}" "$CONTAM_NORMAL_BAM" -o "$POOL.tmp"
  mv "$POOL.tmp" "$POOL"; "$SAMTOOLS" index -@ "$THREADS" "$POOL"
fi
POOL_COV=$("$SAMTOOLS" idxstats "$POOL" | awk -v rl="$READLEN" -v chrs="$CHROM_LIST" \
  'BEGIN{n=split(chrs,a," ");for(i=1;i<=n;i++)keep[a[i]]=1} keep[$1]{len+=$2;mp+=$3} END{if(len>0)printf "%.3f",mp*rl/len; else print 0}')
log "contamination pool mean depth over sim chroms ~= ${POOL_COV}x"

# ----------------------------------------------------------------------------
# 4. mix: impure tumor = p*tumor + (1-p)*normal-pool, at each coverage x purity
# ----------------------------------------------------------------------------
MADE=()
for C in $COVERAGES; do
  for P in $PURITIES; do
    OUT="$RUNDIR/tumor.${C}x.p${P}.bam"
    if [ -f "$OUT" ]; then log "exists -> skip $(basename "$OUT")"; MADE+=("$(basename "$OUT")"); continue; fi
    # sanity: enough normal pool for the (1-p)*C contamination?
    need=$(awk -v p="$P" -v c="$C" 'BEGIN{printf "%.3f",(1-p)*c}')
    awk -v n="$need" -v have="$POOL_COV" 'BEGIN{exit !(n<=have)}' || \
      log "WARN: need ${need}x normal but pool only ${POOL_COV}x -> normal will be capped (purity higher than requested)"
    log "mixing tumor ${C}x purity ${P} -> $(basename "$OUT") ..."
    "$HERE/mix_tumor_normal.sh" -t "$PURE_TUMOR" -n "$POOL" -p "$P" -d "$C" -o "$OUT" \
        --tumor-cov "$TUMOR_SIM_COV" --normal-cov "$POOL_COV" \
        --readlen "$READLEN" --seed "$SEED" --threads "$THREADS"
    MADE+=("$(basename "$OUT")")
  done
done

# ----------------------------------------------------------------------------
# 5. trim intermediates (regeneratable), keep tumor BAMs + truth
# ----------------------------------------------------------------------------
rm -rf "$RUNDIR/fastq_tumor" "$RUNDIR"/bwa_*.log 2>/dev/null || true
log "removed FASTQ/bwa-log intermediates (regeneratable from seed)"

# ----------------------------------------------------------------------------
# 6. README with full provenance
# ----------------------------------------------------------------------------
sv_counts=$(tail -n +2 "$RUNDIR/truth.sv_table.tsv" | cut -f2 | sort | uniq -c | awk '{printf "%s=%s ",$2,$1}')
ntruth=$(($(wc -l < "$RUNDIR/truth.bedpe") - 1))
# tool versions (these CLIs exit non-zero with no args; guard under `set -e`)
artv=$(art_illumina 2>&1 | awk '/Version/{print $0; exit}' || true)
bwav=$(bwa 2>&1 | awk '/Version/{print $2; exit}' || true)
samv=$($SAMTOOLS --version 2>/dev/null | head -1 || true)
cat > "$RUNDIR/README.md" <<EOF
# Simulated tumor panel — $(basename "$RUNDIR")

Generated $(date) on $(hostname) by \`sim/build_sim_panel.sh\` (svaba repo).
**Reproducible: SEED = ${SEED}.**

## What this is
A diploid donor genome carrying **${ntruth} non-overlapping somatic SVs** on
**${CHROMS}** (${sv_counts}), simulated to reads with ART, aligned with BWA-MEM,
then admixed with a **real whole-genome normal** to emulate tumor impurity.
Every simulated SV is somatic (germline-frac 0).

## Tumor BAMs (indexed)
$(for b in "${MADE[@]}"; do echo "- \`$b\`"; done)

Naming: \`tumor.<coverage>x.p<purity>.bam\` — total depth (over sim chroms) = coverage;
\`purity\` = tumor fraction (the rest is real-normal contamination). Also present:
\`tumor.${TUMOR_SIM_COV}x.bam\` = the pure (100%) simulated tumor (intermediate).

## Ground truth (reference coordinates, 0-based half-open)
- \`truth.bedpe\`  — one row per SV: chrom1,start1,end1,chrom2,start2,end2,name,svtype,svlen,somatic,gt,haps
- \`truth.vcf\`    — symbolic <DEL>/<DUP>/<INV>/<INS> + BND records, SOMATIC flag
- \`truth.sv_table.tsv\` — full per-SV dump

## How it was made
1. SVs + donor genome + truth:
   \`\`\`
   python3 simulate_sv_genome.py --ref <REF> --out-dir <here> \\
     --chroms ${CHROMS} --seed ${SEED} ${GEN_ARGS}
   \`\`\`
2. Pure tumor reads: ART \`-ss ${ART_PROFILE} -l ${READLEN} -m ${FRAG_MEAN} -s ${FRAG_SD}\`
   on each haplotype at ${TUMOR_SIM_COV}/2×; BWA-MEM to <REF>; sort/index.
3. Contamination pool: random ${NORMAL_DOWNSAMPLE} downsample (seed ${SEED}) of the
   real normal \`${CONTAM_NORMAL_BAM}\` (mean depth over sim chroms ~= ${POOL_COV}x),
   cached at \`${POOL}\`.
4. Impure tumor = purity·(pure tumor) + (1−purity)·(normal pool), coverage-matched
   over the simulated chromosomes: keep \`p·C/${TUMOR_SIM_COV}\` of tumor +
   \`(1−p)·C/${POOL_COV}\` of pool, force one tumor read-group, merge/sort/index
   (sim/mix_tumor_normal.sh).

## Coordinate / model notes
- The simulated tumor reads are reference-derived (no germline SNPs) **plus** the
  somatic SVs; the real-normal admixture supplies genuine germline background.
  So a het somatic SV (VAF 0.5 in pure tumor) appears at ~0.5·purity in the mix.
- ART emits no PCR/optical duplicates.

## Tools
- ART: ${artv:-?}   |   bwa: ${bwav:-?}   |   ${samv:-samtools ?}

## Run svaba + score
Pair the tumor with your OTHER normal (matched normal for this individual):
\`\`\`
svaba run -t ${RUNDIR}/tumor.${maxC}x.p$(printf '%s' "$PURITIES" | awk '{print $1}').bam -n ${MATCHED_NORMAL_BAM} \\
          -G ${REF} -a eval -p ${THREADS}
svaba postprocess -i eval
# score / visualize: docs/sim_benchmark.html  (load truth.bedpe + eval.bps.sorted.dedup.txt.gz)
#  or: python3 sim/benchmark.py --truth ${RUNDIR}/truth.bedpe --calls eval.bps.sorted.dedup.txt.gz --pass-only --somatic-only
\`\`\`

## Regenerate exactly
\`\`\`
SEED=${SEED} TAG=${TAG} CHROMS="${CHROMS}" COVERAGES="${COVERAGES}" PURITIES="${PURITIES}" \\
NORMAL_DOWNSAMPLE=${NORMAL_DOWNSAMPLE} GEN_ARGS="${GEN_ARGS}" \\
CONTAM_NORMAL_BAM="${CONTAM_NORMAL_BAM}" sim/build_sim_panel.sh
\`\`\`
EOF

log "DONE -> $RUNDIR"
log "  tumor BAMs: ${MADE[*]}"
log "  truth: truth.bedpe / truth.vcf / truth.sv_table.tsv"
log "  see $RUNDIR/README.md"
