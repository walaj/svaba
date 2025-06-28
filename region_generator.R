library(data.table)
options(scipen = 999)

# Load hg38 chromosome sizes from UCSC
chr_sizes <- fread("https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes", col.names = c("chrom", "size"))
chr_sizes <- chr_sizes[chrom %in% paste0("chr", c(1:22, "X"))]

# Centromere positions (approximate from UCSC Table Browser / literature)
centromeres <- fread("
chrom start end
chr1 121700000 124600000
chr2 91800000 94000000
chr3 87800000 94000000
chr4 48200000 50700000
chr5 46100000 50700000
chr6 58600000 63400000
chr7 58000000 62700000
chr8 43900000 48100000
chr9 47100000 50700000
chr10 39100000 42900000
chr11 51400000 55800000
chr12 34700000 38600000
chr13 16000000 19500000
chr14 16000000 19000000
chr15 17000000 20500000
chr16 35300000 38200000
chr17 22100000 25800000
chr18 15400000 18700000
chr19 24400000 28000000
chr20 25800000 30000000
chr21 10900000 14300000
chr22 12000000 15400000
chrX 58000000 61900000
")

# Construct p and q arms
bed <- rbindlist(lapply(1:nrow(chr_sizes), function(i) {
  chrom <- chr_sizes$chrom[i]
  size  <- chr_sizes$size[i]
  cent  <- centromeres[chrom == chrom]

  if (nrow(cent) == 0) return(NULL)

  data.table(
    chrom = c(chrom, chrom),
    start = c(0, cent$end),
    end   = c(cent$start, size)
  )
}))

# Natural chromosome order: chr1 chr22, chrX
bed[, chrom := factor(chrom, levels = paste0("chr", c(1:22, "X")))]

# Sort by chromosome and start position
setorder(bed, chrom, start)

# Write to BED file
fwrite(bed, "hg38_arms_excl_centromeres.bed", sep = "\t", col.names = FALSE)
