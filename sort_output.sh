#!/bin/bash

SAM=samtools
ID=$1
MEM=16G
THREADS=${2:-4}

# Suffixes that need deduplication (window overlap can produce exact duplicate reads)
DEDUP_SUFFIXES="weird corrected discordant"

# If set to "1", run an extra module that splits the discordant/weird/corrected
# BAMs into per-source BAMs, using the first 4 chars of each read name as the
# source tag (e.g. a read named "t001_H1234..." gets routed into
# ${ID}.${suffix}.t001.bam). Override at invocation with:
#     SPLIT_BY_SOURCE=1 ./sort_output.sh ID
SPLIT_BY_SOURCE="${SPLIT_BY_SOURCE:-1}"

# Preprocessing: merge thread BAMs if necessary
for suffix in discordant weird corrected; do
  pattern="${ID}.thread*.${suffix}.bam"
  target="${ID}.${suffix}.bam"
  shopt -s nullglob
  bam_files=($pattern)
  shopt -u nullglob
  if [[ ${#bam_files[@]} -gt 1 ]]; then
    echo "Merging ${#bam_files[@]} BAM files for suffix '$suffix' into $target"
    ${SAM} merge -f -@ ${THREADS} "$target" "${bam_files[@]}" && rm "${bam_files[@]}"
  elif [[ ${#bam_files[@]} -eq 1 ]]; then
    echo "Renaming single-thread BAM ${bam_files[0]} to $target"
    mv "${bam_files[0]}" "$target"
  fi
done

# Now process all expected BAMs for sorting, dedup, and indexing
for suffix in weird corrected contigs discordant; do
  bam="${ID}.${suffix}.bam"
  sorted="${ID}.${suffix}.sorted.bam"

  if [[ ! -f "$bam" ]]; then
    echo "Skipping: $bam not found."
    continue
  fi

  echo "Sorting $bam..."
  ${SAM} sort -@ ${THREADS} -m ${MEM} -o "$sorted" "$bam"
  mv "$sorted" "$bam"

  # Deduplicate by read name + flags for suffixes with window overlap
  if [[ " ${DEDUP_SUFFIXES} " == *" ${suffix} "* ]]; then
    echo "Deduplicating $bam (removing exact read name + flag duplicates)..."
    deduped="${ID}.${suffix}.deduped.bam"
    ${SAM} view -@ ${THREADS} -h "$bam" \
      | awk 'BEGIN{OFS="\t"} /^@/{print; next} {key=$1"\t"$2; if(!seen[key]++){print}}' \
      | ${SAM} view -@ ${THREADS} -bS -o "$deduped" -
    mv "$deduped" "$bam"
  fi

  ${SAM} index -@ ${THREADS} "$bam"
  echo "Done: $bam"
done

# Optional: split BAMs by read-name source prefix (first 4 chars of QNAME).
# Enabled by setting SPLIT_BY_SOURCE=1 at the top of this script or on the
# command line. Produces files like ${ID}.${suffix}.${prefix}.bam (e.g.
# ${ID}.weird.t001.bam) for each distinct 4-char prefix found.
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

    # Reassemble each per-prefix SAM into a sorted+indexed BAM.
    shopt -s nullglob
    for body in "${tmpdir}"/body.*.sam; do
      base=$(basename "$body" .sam)
      prefix="${base#body.}"
      out="${ID}.${suffix}.${prefix}.bam"
      echo "    Writing $out"
      cat "$header" "$body" | ${SAM} view -bS - | \
        ${SAM} sort -m ${MEM} -o "$out" - && \
        ${SAM} index "$out"
    done
    shopt -u nullglob

    rm -rf "$tmpdir"
  done
fi
