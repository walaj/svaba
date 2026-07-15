#!/usr/bin/env bash
# sim_config.sh -- central configuration for the svaba SV simulator.
# Sourced by run_sim.sh and mix_tumor_normal.sh. Override any of these by
# exporting them before invoking, or edit here. All paths may be relative.

# ---- reference + output ----------------------------------------------------
: "${REF:=$HOME/ref_genome/Homo_sapiens_assembly38.fasta}"   # must be bwa-indexed + .fai
: "${OUTDIR:=$PWD/sim_out}"                                   # all outputs land here
: "${SEED:=42}"
: "${THREADS:=8}"

# Which chromosomes to mutate + simulate from. Reads still align to the FULL
# reference. Use a subset (e.g. "chr20-22") for a fast smoke test.
: "${CHROMS:=chr1-22,chrX,chrY}"

# Optional BED of regions to AVOID placing SVs in (centromeres / blacklist).
# Recommended: the svaba combined blacklist.
: "${AVOID_BED:=}"

# ---- ART read simulation ---------------------------------------------------
: "${ART_PROFILE:=HS25}"     # HiSeq 2500 150bp error profile (art_illumina -ss)
: "${READLEN:=150}"
: "${FRAG_MEAN:=400}"        # mean fragment (insert) size
: "${FRAG_SD:=50}"

# Maximum tumor coverage to SIMULATE (once); lower coverages are produced by
# downsampling the aligned BAM (fast, no re-simulation).
: "${MAXCOV:=30}"
: "${COVERAGES:=5 15 30}"    # target tumor coverages to emit (each <= MAXCOV)

# ---- normal sample ---------------------------------------------------------
# For the tumor/normal pair and for impurity admixture you need a normal BAM.
#   * Default: use a REAL blood normal (recommended; all simulated SVs are then
#     unambiguously somatic against it).
#   * Or set SIM_NORMAL=1 to simulate a clean normal from the reference instead.
: "${NORMAL_BAM:=$HOME/Downloads/svaba_compare/blood.recal.bam}"
: "${SIM_NORMAL:=0}"         # 1 => simulate a clean normal from REF (ignores NORMAL_BAM)
: "${NORMAL_COV:=30}"        # coverage to simulate the normal at (if SIM_NORMAL=1)

# ---- SV generator knobs (passed to simulate_sv_genome.py) ------------------
# Counts/sizes that yield ~3,200 well-spaced SVs across the genome. Scale up
# --n-* for denser coverage. --min-gap keeps each event in its own ~independent
# locus (>= svaba's assembly window). Default is all-somatic (--germline-frac 0).
: "${GEN_ARGS:=--n-del-small 800 --n-ins-small 800 --n-del 600 --n-dup 400 --n-inv 400 --n-tra 200 --min-gap 50000 --hom-frac 0.3 --germline-frac 0.0}"

# ---- tool binaries (override if not on PATH) -------------------------------
: "${ART:=art_illumina}"
: "${BWA:=bwa}"
: "${SAMTOOLS:=samtools}"
: "${PYTHON:=python3}"
: "${SIM_PY:=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/simulate_sv_genome.py}"
