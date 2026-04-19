#!/usr/bin/env bash
#
# combine_blacklists.sh — combine multiple BED blacklists into one file.
#
# "Combining" blacklists means a UNION over intervals: a region ends up in the
# output if ANY of the input files covered it. This is what you almost always
# want when you're assembling a single --blacklist argument for svaba run out
# of several sources (high-signal regions, low-mappability, simple repeats,
# non-standard contigs, high-runtime regions, ad-hoc bad list, etc.).
#
# (If you actually need the set intersection — regions blacklisted by ALL
# sources simultaneously — this is not that script. For that, use
# `bedtools multiinter -i FILES | awk '$4==NFILES'` or iterated
# `bedtools intersect`. Intersect-of-blacklists is almost never the right
# semantics; you'd end up with a tiny output.)
#
# What this script does:
#   1. cat every input BED to stdout with a 5th-column source tag = basename.
#      Any existing 4th-column label is preserved.
#   2. Filter out obvious non-BED noise: empty lines, `#` comments,
#      `track`/`browser` UCSC header lines.
#   3. By default: sort by (chrom, start) and emit as a 5-column BED.
#   4. With --merge: run `bedtools merge` to collapse overlapping/adjacent
#      intervals into one, and collapse the set of contributing source tags
#      into a comma-separated list in column 4.
#
# Usage:
#     ./combine_blacklists.sh -o OUT.bed [--merge [--slop N]] FILE1.bed FILE2.bed ...
#
# Flags:
#     -o, --output OUT     output path (required)
#     -m, --merge          sort + bedtools merge the concatenation (union of
#                          overlapping intervals, with label-distinct column).
#                          Requires `bedtools` on PATH.
#     -d, --merge-dist N   pass to bedtools merge -d (intervals within N bp
#                          are merged). Default 0 (touching/overlapping only).
#         --slop N         grow each interval by N bp on both sides before
#                          merging (expands the blacklist). Default 0.
#                          Requires `bedtools` and --merge.
#     -c, --clip GENOME    clip every interval's end to the real contig length
#                          in GENOME, and drop intervals whose start is past the
#                          contig. GENOME may be a .fai, a two-column tsv
#                          (chrom<TAB>length), or svaba's tracks/chr fasta-
#                          header dump (>chrN AC:...  LN:... ...). Unknown
#                          contigs are passed through unchanged (with a
#                          warning). This is what saves you from a 250,000,000-
#                          end sentinel in one of the input BEDs blowing the
#                          reported covered-bp up to absurd levels.
#         --no-label       drop the label column(s); emit 3-col BED only
#     -q, --quiet          suppress progress to stderr
#     -h, --help           this message
#
# Typical invocation (regenerate the project's combined blacklist):
#
#   ./combine_blacklists.sh --merge \
#     -o tracks/hg38.combined_blacklist.bed \
#     tracks/hg38.blacklist.bed \
#     tracks/hg38.high_runtime.bed \
#     tracks/hg38.extra_badlist.bed \
#     tracks/hg38.nonstd_chr.blacklist.bed \
#     tracks/hg38.rmsk.simple_repeat.bed
#
# The script is idempotent — re-running it always rebuilds OUT from current
# inputs. Change `high_runtime.bed`, rerun, and the combined file updates.

set -euo pipefail

OUT=""
MERGE=0
MERGE_DIST=0
SLOP=0
CLIP=""
NO_LABEL=0
QUIET=0
INPUTS=()

log() { [[ $QUIET -eq 1 ]] || echo "$@" >&2; }

usage() { sed -n '2,60p' "$0" | sed 's/^# \{0,1\}//'; exit "${1:-0}"; }

# ---- arg parsing ----
while (( $# )); do
  case "$1" in
    -o|--output)      OUT="$2"; shift 2 ;;
    -m|--merge)       MERGE=1; shift ;;
    -d|--merge-dist)  MERGE_DIST="$2"; shift 2 ;;
    --slop)           SLOP="$2"; shift 2 ;;
    -c|--clip)        CLIP="$2"; shift 2 ;;
    --no-label)       NO_LABEL=1; shift ;;
    -q|--quiet)       QUIET=1; shift ;;
    -h|--help)        usage 0 ;;
    --)               shift; while (( $# )); do INPUTS+=("$1"); shift; done ;;
    -*)               echo "Unknown flag: $1" >&2; usage 2 ;;
    *)                INPUTS+=("$1"); shift ;;
  esac
done

[[ -z "$OUT" ]]       && { echo "ERROR: -o OUTPUT is required" >&2; usage 2; }
[[ ${#INPUTS[@]} -lt 1 ]] && { echo "ERROR: no input BED files given" >&2; usage 2; }

# ---- sanity on inputs ----
for f in "${INPUTS[@]}"; do
  [[ -f "$f" ]] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done

if [[ $MERGE -eq 1 || $SLOP -gt 0 ]]; then
  command -v bedtools >/dev/null 2>&1 \
    || { echo "ERROR: --merge/--slop require 'bedtools' on PATH" >&2; exit 1; }
fi

if [[ $SLOP -gt 0 && $MERGE -ne 1 ]]; then
  echo "ERROR: --slop requires --merge (expanded intervals need to be merged to stay collapsed)" >&2
  exit 2
fi

# ---- work in a tempdir ----
TMPDIR_LOCAL="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT
CAT="$TMPDIR_LOCAL/cat.bed"
SORTED="$TMPDIR_LOCAL/sorted.bed"
GENOME_TSV="$TMPDIR_LOCAL/genome.tsv"

# ---- optional: parse clip-genome file into a flat chrom<TAB>length TSV ----
#
# Accepts three input shapes; we auto-detect:
#   .fai       : name<TAB>len<TAB>offset<TAB>linebases<TAB>linewidth
#   genome tsv : name<TAB>len                                (optional extra cols ignored)
#   tracks/chr : >chrN  AC:...  LN:NNNN  ...                 (fasta-header dump)
if [[ -n "$CLIP" ]]; then
  [[ -f "$CLIP" ]] || { echo "ERROR: --clip file not found: $CLIP" >&2; exit 1; }
  awk -v OFS='\t' '
    /^>/ {
      name=$1; sub(/^>/, "", name)
      ln=""
      for (i=2;i<=NF;i++) if ($i ~ /^LN:/) { ln=substr($i,4); break }
      if (ln!="" && ln+0>0) print name, ln
      next
    }
    /^[[:space:]]*$/ {next}
    /^[[:space:]]*#/ {next}
    {
      # .fai or 2-col genome. col1=name, col2=length.
      if ($2 ~ /^[0-9]+$/ && $2+0 > 0) print $1, $2
    }
  ' "$CLIP" > "$GENOME_TSV"
  CLIP_N="$(wc -l < "$GENOME_TSV" | awk '{print $1}')"
  [[ "$CLIP_N" -gt 0 ]] || { echo "ERROR: --clip file $CLIP has no usable (name, length) entries" >&2; exit 1; }
  log "combine_blacklists.sh: clip genome loaded from $CLIP ($CLIP_N contigs)"
fi

log "combine_blacklists.sh: building $OUT from ${#INPUTS[@]} source(s)"

# ---- step 1: concatenate with a 5th-column source tag ----
: > "$CAT"
for f in "${INPUTS[@]}"; do
  src="$(basename "$f")"
  # Drop blank / comment / UCSC-header lines. For each record, collapse leading
  # whitespace and ensure we end up tab-separated. Keep cols 1-3, promote any
  # existing col 4 to the label, and append src as col 5.
  awk -v src="$src" '
    BEGIN{OFS="\t"}
    /^[[:space:]]*$/ {next}
    /^[[:space:]]*#/ {next}
    /^[[:space:]]*(track|browser)[[:space:]]/ {next}
    {
      # Normalize: the input may already be tab-separated (BED spec) or
      # whitespace-separated. Use awk default FS (runs of whitespace) but
      # write tabs on output.
      chr=$1; start=$2; end=$3
      if (start !~ /^[0-9]+$/ || end !~ /^[0-9]+$/) next
      label=""
      if (NF >= 4) {
        label=$4
        for (i=5;i<=NF;i++) label=label" "$i
      }
      if (label=="") print chr, start, end, "-", src
      else          print chr, start, end, label, src
    }
  ' "$f" >> "$CAT"
  log "  + $(printf "%-45s" "$src")  $(wc -l < "$f" | awk '{print $1}') lines"
done

RAW_LINES="$(wc -l < "$CAT" | awk '{print $1}')"
log "  = $RAW_LINES raw intervals before sort/merge"

# ---- step 2: sort ----
# LC_ALL=C to keep sort output stable across locales and compatible with
# bedtools' lexicographic chrom ordering.
LC_ALL=C sort -k1,1 -k2,2n "$CAT" > "$SORTED"

# ---- step 2b: optional clip-to-genome ----
#
# Applied BEFORE merge and slop so downstream steps see real coordinates.
# Intervals on unknown contigs are passed through unchanged (with a single
# warning). Intervals whose start is beyond the contig length are dropped;
# intervals whose end overshoots are clipped to the contig length.
if [[ -n "$CLIP" ]]; then
  CLIPPED="$TMPDIR_LOCAL/clipped.bed"
  awk -v OFS='\t' -v genomefile="$GENOME_TSV" '
    BEGIN {
      while ((getline < genomefile) > 0) { len[$1] = $2 }
      close(genomefile)
      dropped=0; clipped=0; passthrough=0; unknown_n=0
    }
    {
      chr=$1; s=$2; e=$3
      if (!(chr in len)) {
        if (!(chr in unknown_seen)) { unknown_seen[chr]=1; unknown_n++ }
        passthrough++
        print
        next
      }
      L = len[chr]
      if (s+0 >= L) { dropped++; next }
      if (e+0 > L)  { e = L; clipped++ }
      $2=s; $3=e
      print
    }
    END {
      fmt="  · clip: %d clipped to contig length, %d dropped past contig end, %d intervals on %d unknown contig(s) passed through\n"
      printf fmt, clipped, dropped, passthrough, unknown_n > "/dev/stderr"
      if (unknown_n > 0) {
        n=0
        for (c in unknown_seen) {
          if (n < 5) printf "    · unknown contig example: %s\n", c > "/dev/stderr"
          n++
        }
        if (unknown_n > 5) printf "    · ... (%d unknown contigs total)\n", unknown_n > "/dev/stderr"
      }
    }
  ' "$SORTED" > "$CLIPPED"
  SORTED="$CLIPPED"
fi

# ---- step 3: emit ----
if [[ $MERGE -eq 1 ]]; then
  # bedtools merge collapses overlapping / touching intervals. We carry the
  # source tag (col 5) into the collapsed output as a comma-separated set
  # of contributing sources via -c 5 -o distinct; and also collapse the
  # original labels via col 4 so downstream tools still see a human-readable
  # label column.
  #
  # Layout of the intermediate we feed to bedtools:
  #   chr start end label src
  # After merge with -c 4,5 -o distinct,distinct we get:
  #   chr start end label_set src_set

  if [[ $SLOP -gt 0 ]]; then
    GEN="$TMPDIR_LOCAL/genome.tsv"
    # Try to derive a genome file from the first input's chroms as max(end).
    # This is a permissive fallback; if you have a real .fai / .genome file
    # available, you can swap this with `bedtools slop -g your.genome`.
    awk 'BEGIN{OFS="\t"} {if($3>m[$1]) m[$1]=$3} END{for (k in m) print k, m[k]}' \
      "$SORTED" > "$GEN"
    SLOPPED="$TMPDIR_LOCAL/slopped.bed"
    bedtools slop -i "$SORTED" -g "$GEN" -b "$SLOP" > "$SLOPPED"
    LC_ALL=C sort -k1,1 -k2,2n "$SLOPPED" > "$SORTED"
    log "  · slopped +/-${SLOP}bp (genome inferred from max end per chrom)"
  fi

  MERGED="$TMPDIR_LOCAL/merged.bed"
  bedtools merge -i "$SORTED" -d "$MERGE_DIST" -c 4,5 -o distinct,distinct \
    > "$MERGED"
  MERGED_LINES="$(wc -l < "$MERGED" | awk '{print $1}')"
  log "  · bedtools merge (d=$MERGE_DIST) : $RAW_LINES -> $MERGED_LINES intervals"

  if [[ $NO_LABEL -eq 1 ]]; then
    cut -f1-3 "$MERGED" > "$OUT"
  else
    # 5 cols: chr start end label-set source-set
    cp "$MERGED" "$OUT"
  fi
else
  if [[ $NO_LABEL -eq 1 ]]; then
    cut -f1-3 "$SORTED" > "$OUT"
  else
    # 5 cols: chr start end label src (pre-merge, duplicates preserved)
    cp "$SORTED" "$OUT"
  fi
fi

FINAL_LINES="$(wc -l < "$OUT" | awk '{print $1}')"
FINAL_BP="$(awk 'BEGIN{s=0} {s+=$3-$2} END{print s+0}' "$OUT")"
log "combine_blacklists.sh: wrote $OUT  ($FINAL_LINES intervals, $(printf "%'d" "$FINAL_BP") bp covered)"

# Sanity check: the human genome is ~3.1 Gbp (~3.2 Gbp incl. alts). If the
# reported coverage is larger than a loose 10 Gbp threshold, the input almost
# certainly contains synthetic "end-of-contig" sentinels (e.g. end=250000000
# on short alt contigs) and the user should rerun with --clip GENOME to get
# a meaningful number.
if [[ $FINAL_BP -gt 10000000000 ]]; then
  log ""
  log "WARNING: $(printf "%'d" $FINAL_BP) bp is larger than any reasonable reference."
  if [[ -z "$CLIP" ]]; then
    log "  Likely cause: an input BED uses oversized end coordinates as a 'whole"
    log "  contig' shorthand (e.g. end=250000000 on short alt contigs). Rerun"
    log "  with  --clip tracks/chr  (or a .fai / genome tsv) to clip each"
    log "  interval to its real contig length before the bp tally."
  else
    log "  --clip was used; the oversize is probably legitimate (very large"
    log "  genome) or caused by intervals on contigs missing from your clip"
    log "  genome file — check the 'unknown contig' lines above."
  fi
fi
