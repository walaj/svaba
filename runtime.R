#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
})

# Parse arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: ./runtime_plot.R input.tsv output.pdf")
}

input_file  <- args[1]
output_file <- args[2]

# Read data
dt <- fread(input_file)

# Only keep chromosomes 1-22 and X
allowed_chr <- paste0("chr", c(1:22, "X"))
dt <- dt[chromosome %in% allowed_chr]

# Convert start to Mb
dt[, start_mb := start / 1e6]

# Region label for top runtime labels
dt[, region_label := paste0(chromosome, ":", start, "-", end)]

# For each chromosome, find top 10 regions by runtime
top10 <- dt[, .SD[order(-runtime_seconds)][1:10], by = chromosome]

# Plot with facet by chromosome
p <- ggplot(dt, aes(x = start_mb, y = runtime_seconds)) +
  geom_point(size = 0.7) +
  geom_line(aes(group = 1), alpha = 0.7) +  # explicitly group the line
  geom_text_repel(
    data         = top10,
    aes(label = region_label),
    size         = 2,
    max.overlaps = Inf,
    box.padding  = 0.5
  ) +
  facet_wrap(~chromosome, scales = "free_x", ncol = 5) +
  labs(
    title = "Runtime per Genomic Region by Chromosome",
    x     = "Genomic Start Position (Mb)",
    y     = "Wall Time (seconds)"
  ) +
  theme_minimal(base_size = 9) +
  theme(strip.text = element_text(face = "bold")) +
  scale_x_continuous(expand = expansion(add = c(0.01, 0.01)))  # avoid zero-width panels

# Save a large PDF
ggsave(output_file, plot = p, width = 20, height = 16, dpi = 300)
