#!/usr/bin/env bash
#
# filter_contig_supporting_reads.sh
#
# Filter a svaba *.corrected.bam down to just the reads whose r2c (read-to-
# contig) BWA alignment landed on a single named contig. The contig name is
# matched against the `bi` aux tag that svaba writes onto each corrected
# read (one tag value per r2c hit, comma-separated). Output is a coordinate-
# sorted+indexed BAM for easy IGV inspection.

set -euo pipefail

PROG="$(basename "$0")"

usage() {
  cat <<EOF
Usage: ${PROG} -i ID -c CONTIG [-o OUT] [-t THREADS]
       ${PROG} ID CONTIG                       # positional shorthand

Filter \${ID}.corrected.bam down to reads that r2c-aligned to CONTIG (matched
against the 'bi' aux tag) and write a sorted+indexed BAM.

Required:
  -i, --id      ID         svaba sample/run prefix; expects \${ID}.corrected.bam
                           in the current directory (or a path; .corrected.bam
                           is appended if it doesn't end in .bam).
  -c, --contig  CONTIG     Contig name to filter on (e.g.
                           c_fermi_chr12_13405500_13430500_13C).

Optional:
  -o, --out     OUT        Output BAM path (default: \${ID}.\${CONTIG}.bam).
  -t, --threads N          Threads passed to samtools sort (default: 4).
  -h, --help               Show this help and exit.

Example:
  ${PROG} -i tumor_run1 -c c_fermi_chr12_13405500_13430500_13C
  ${PROG} tumor_run1 c_fermi_chr12_13405500_13430500_13C -o chr12_13C.bam
EOF
}

# --- arg parsing ------------------------------------------------------------
ID=""
CONTIG=""
OUT=""
THREADS=4
POS=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    -h|--help)    usage; exit 0 ;;
    -i|--id)      ID="$2";      shift 2 ;;
    -c|--contig)  CONTIG="$2";  shift 2 ;;
    -o|--out)     OUT="$2";     shift 2 ;;
    -t|--threads) THREADS="$2"; shift 2 ;;
    --)           shift; POS+=("$@"); break ;;
    -*)           echo "${PROG}: unknown option: $1" >&2; usage >&2; exit 2 ;;
    *)            POS+=("$1"); shift ;;
  esac
done

# Positional fallback: ID CONTIG
if [[ -z "$ID"     && ${#POS[@]} -ge 1 ]]; then ID="${POS[0]}"; fi
if [[ -z "$CONTIG" && ${#POS[@]} -ge 2 ]]; then CONTIG="${POS[1]}"; fi

if [[ -z "$ID" || -z "$CONTIG" ]]; then
  echo "${PROG}: missing required argument(s)" >&2
  usage >&2
  exit 2
fi

# --- resolve input/output paths --------------------------------------------
if [[ "$ID" == *.bam ]]; then
  IN_BAM="$ID"
  ID_BASE="${ID%.corrected.bam}"
  ID_BASE="${ID_BASE%.bam}"
else
  IN_BAM="${ID}.corrected.bam"
  ID_BASE="$ID"
fi

if [[ ! -f "$IN_BAM" ]]; then
  echo "${PROG}: input BAM not found: ${IN_BAM}" >&2
  exit 1
fi

if [[ -z "$OUT" ]]; then
  OUT="${ID_BASE}.${CONTIG}.bam"
fi

# --- check tools ------------------------------------------------------------
command -v samtools >/dev/null 2>&1 || {
  echo "${PROG}: samtools not found on PATH" >&2; exit 127;
}

# --- filter -----------------------------------------------------------------
echo "${PROG}: ID        = ${ID_BASE}" >&2
echo "${PROG}: input     = ${IN_BAM}"  >&2
echo "${PROG}: contig    = ${CONTIG}"  >&2
echo "${PROG}: output    = ${OUT}"     >&2
echo "${PROG}: threads   = ${THREADS}" >&2

samtools view -h -e "[bi]=~\"${CONTIG}\"" "${IN_BAM}" \
  | samtools sort -@"${THREADS}" -o "${OUT}" -
samtools index "${OUT}"

echo "${PROG}: wrote ${OUT} (+ ${OUT}.bai)" >&2
