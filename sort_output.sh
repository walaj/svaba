#!/bin/bash

SAM=/usr/bin/samtools
ID=$1
MEM=16G

# Preprocessing: merge thread BAMs if necessary
for suffix in weird corrected; do
  pattern="${ID}.thread*.${suffix}.bam"
  target="${ID}.${suffix}.bam"

  shopt -s nullglob
  bam_files=($pattern)
  shopt -u nullglob

  if [[ ${#bam_files[@]} -gt 1 ]]; then
      echo "Merging ${#bam_files[@]} BAM files for suffix '$suffix' into $target"
      echo "${SAM} merge -f $target ${bam_files[@]}"
    ${SAM} merge -f "$target" "${bam_files[@]}"
  elif [[ ${#bam_files[@]} -eq 1 ]]; then
    echo "Renaming single-thread BAM ${bam_files[0]} to $target"
    mv "${bam_files[0]}" "$target"
  fi
done

# Now process all expected BAMs for sorting and indexing
for suffix in weird corrected contigs; do
  bam="${ID}.${suffix}.bam"
  sorted="${ID}.${suffix}.sorted.bam"

  echo "Processing $bam..."

  if [[ ! -f "$bam" ]]; then
    echo "  Skipping: $bam not found."
    continue
  fi

  ${SAM} sort -m ${MEM} -o "$sorted" "$bam" && \
    mv "$sorted" "$bam" && \
    ${SAM} index "$bam"

  echo "  Done: $bam"
done
