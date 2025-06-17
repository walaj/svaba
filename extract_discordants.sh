#!/bin/bash

prefix="${1%%.*}"
/usr/bin/samtools view -h $1 | grep -E '^@|DC:Z' | /usr/bin/samtools view -b -o tmp.filtered.bam
/usr/bin/samtools sort -o ${prefix}.dc.bam tmp.filtered.bam
/usr/bin/samtools index ${prefix}.dc.bam
rm tmp.filtered.bam
