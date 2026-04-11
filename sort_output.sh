#!/bin/bash

SAM=samtools
ID=$1
MEM=16G
THREADS=${2:-4}

# Suffixes that need deduplication (window overlap can produce exact duplicate reads)
DEDUP_SUFFIXES="weird corrected discordant"

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
