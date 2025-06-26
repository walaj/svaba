#!/bin/bash

prefix=$1
##prefix="${1%%.*}"
/usr/bin/samtools view -h ${prefix}.weird.bam | grep -E '^@|DC:Z' | /usr/bin/samtools view -b -o tmp.filtered.bam
/usr/bin/samtools sort -o ${prefix}.dc.bam tmp.filtered.bam
/usr/bin/samtools index ${prefix}.dc.bam
rm tmp.filtered.bam

TARGET=$2
if [[ -n "$TARGET" ]]; then
    /usr/bin/samtools view -h ${prefix}.dc.bam | \
	awk -v pat="${TARGET}" 'BEGIN{OFS="\t"} /^@/ {print; next} $0 ~ pat {print}' | \
	/usr/bin/samtools view -b -o dc.${TARGET}.bam
  /usr/bin/samtools index dc.${TARGET}.bam
fi
