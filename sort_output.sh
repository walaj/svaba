#!/bin/bash
#
# sort_output.sh — post-process svaba's per-suffix output BAMs.
#
# The heavy lifting (coord-sort + streaming dedup of exact (qname, flag)
# duplicates from overlapping assembly windows) now lives in the C++
# `svaba postprocess` subcommand, which avoids the old
# `samtools view | awk | samtools view` round-trip through SAM text and
# keeps the dedup hash locus-local rather than global.
#
# This script is the outer glue that:
#   1. Merges per-thread BAMs (${ID}.thread*.${suffix}.bam) into a single
#      ${ID}.${suffix}.bam per suffix, if needed.
#   2. Delegates sort + dedup to `svaba postprocess`.
#   3. Indexes the final BAMs.
#   4. Optionally splits by read-name source prefix (SPLIT_BY_SOURCE=1).
#
# Usage:
#     ./sort_output.sh <ID> [THREADS]
#     SPLIT_BY_SOURCE=1 ./sort_output.sh <ID> [THREADS]

set -e

SAM=samtools
SVABA=${SVABA:-~/git/svaba/build/svaba}       # override if not on PATH
ID=$1
THREADS=${2:-4}
MEM=${MEM:-2G}              # per-samtools-sort-thread memory

if [[ -z "$ID" ]]; then
  echo "Usage: $0 <ID> [THREADS]" >&2
  exit 2
fi

# ---------------------------------------------------------------------------
# Step 1: Merge thread BAMs.
#
# If svaba was run with multiple threads, each suffix may be spread across
# ${ID}.thread0.${suffix}.bam, ${ID}.thread1.${suffix}.bam, etc. Coalesce
# them into a single ${ID}.${suffix}.bam.
# ---------------------------------------------------------------------------
for suffix in discordant weird corrected; do
  target="${ID}.${suffix}.bam"
  shopt -s nullglob
  bam_files=("${ID}".thread*."${suffix}".bam)
  shopt -u nullglob

  if [[ ${#bam_files[@]} -gt 1 ]]; then
    echo "Merging ${#bam_files[@]} thread BAMs for '$suffix' -> $target"
    ${SAM} merge -f -@ "${THREADS}" "$target" "${bam_files[@]}"
    rm -f "${bam_files[@]}"
  elif [[ ${#bam_files[@]} -eq 1 ]]; then
    echo "Renaming single-thread BAM ${bam_files[0]} -> $target"
    mv "${bam_files[0]}" "$target"
  fi
done

# ---------------------------------------------------------------------------
# Step 2: Sort + streaming dedup via `svaba postprocess`.
#
# `svaba postprocess` is a thin C++ orchestrator: for each suffix it shells
# out to `samtools sort` (unchanged — that's where the memory-limited,
# multithreaded external sort still lives) and then runs a native
# SeqLib-based streaming dedup over the sorted BAM. The dedup is the part
# that's actually new in C++: it replaces the old
# `samtools view | awk | samtools view` text round-trip with a single pass
# over BAM records keyed on (qname, flag) with a locus-local hash set.
#
# Handles {weird, corrected, discordant, contigs} in parallel, splits the
# thread budget across concurrent per-suffix jobs (passed through to
# `samtools sort -@` within each job), dedups (qname, flag) duplicates for
# the three dedup-eligible suffixes, leaves `contigs` alone except for
# sorting. Progress prints per 1M reads per suffix during dedup.
# ---------------------------------------------------------------------------
echo "Running '${SVABA} postprocess' (samtools sort + native dedup) with $THREADS threads..."
${SVABA} postprocess -i "$ID" -t "$THREADS" -m "$MEM"

# ---------------------------------------------------------------------------
# Step 3: Index the final BAMs.
#
# Kept here rather than in the C++ module because indexing is cheap, single-
# pass, and not a bottleneck; samtools does it fine.
# ---------------------------------------------------------------------------
for suffix in weird corrected discordant contigs; do
  bam="${ID}.${suffix}.bam"
  [[ -f "$bam" ]] || continue
  echo "Indexing $bam"
  ${SAM} index -@ "${THREADS}" "$bam"
done

# ---------------------------------------------------------------------------
# Step 4 (optional): Split BAMs by read-name source prefix.
#
# Enabled by SPLIT_BY_SOURCE=1. Routes each record into
# ${ID}.${suffix}.${prefix}.bam based on the first 4 chars of QNAME
# (the "source tag" like t001 / n001), then sorts+indexes each output.
#
# This step is still shell because it's a one-shot demux that happens to be
# tolerable via awk+samtools. If it becomes a bottleneck, moving it into
# `svaba postprocess --split-by-source` is a mechanical port: open N
# BamWriters keyed by the first 4 chars, route each record on the fly.
# ---------------------------------------------------------------------------
SPLIT_BY_SOURCE="${SPLIT_BY_SOURCE:-0}"

if [[ "$SPLIT_BY_SOURCE" == "1" ]]; then
  echo "SPLIT_BY_SOURCE=1: splitting BAMs by read-name source prefix"
  for suffix in discordant weird corrected; do
    bam="${ID}.${suffix}.bam"
    if [[ ! -f "$bam" ]]; then
      echo "  Skipping split: $bam not found."
      continue
    fi

    echo "  Splitting $bam by first 4 chars of QNAME..."
    tmpdir=$(mktemp -d)
    header="${tmpdir}/header.sam"
    ${SAM} view -H "$bam" > "$header"

    # Single pass: route each record to a per-prefix body file.
    ${SAM} view "$bam" | awk -v tmpdir="$tmpdir" '{
      prefix = substr($1, 1, 4)
      out = tmpdir "/body." prefix ".sam"
      print >> out
    }'

    # Reassemble each per-prefix SAM into a sorted+indexed BAM, in parallel.
    shopt -s nullglob
    bodies=("${tmpdir}"/body.*.sam)
    n_bodies=${#bodies[@]}
    if [[ $n_bodies -gt 0 ]]; then
      per_job_threads=$(( THREADS / n_bodies ))
      [[ $per_job_threads -lt 1 ]] && per_job_threads=1
      max_parallel=$(( THREADS < n_bodies ? THREADS : n_bodies ))
      echo "    Reassembling $n_bodies per-prefix BAMs ($max_parallel parallel, $per_job_threads samtools threads each)"

      pids=()
      for body in "${bodies[@]}"; do
        while [[ ${#pids[@]} -ge $max_parallel ]]; do
          new_pids=()
          for pid in "${pids[@]}"; do
            kill -0 "$pid" 2>/dev/null && new_pids+=("$pid")
          done
          pids=("${new_pids[@]}")
          [[ ${#pids[@]} -ge $max_parallel ]] && sleep 0.1
        done

        base=$(basename "$body" .sam)
        prefix="${base#body.}"
        out="${ID}.${suffix}.${prefix}.bam"
        echo "      Writing $out"
        (
          cat "$header" "$body" \
            | ${SAM} view -@ ${per_job_threads} -u -bS - \
            | ${SAM} sort -@ ${per_job_threads} -m ${MEM} -o "$out" - \
            && ${SAM} index -@ ${per_job_threads} "$out"
        ) &
        pids+=("$!")
      done
      wait
    fi
    shopt -u nullglob

    rm -rf "$tmpdir"
  done
fi

echo "sort_output.sh: done."
