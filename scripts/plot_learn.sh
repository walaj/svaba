#!/usr/bin/env bash
# plot_learn.sh — plot per-RG insert-size distributions from svaba learn data
#
# Usage:
#   plot_learn.sh <ID_or_file> [output.pdf]
#
# If the argument ends in .tsv or .tsv.gz, treat it as the learn file directly.
# Otherwise treat it as a run ID and glob for ${ID}.*.learn.tsv.gz in the
# current directory.
#
# Output PDF defaults to ${ID}.learn.pdf or ${stem}.learn.pdf.
#
# Requires: Rscript, R packages ggplot2 + data.table (both on CRAN).

set -euo pipefail

usage() {
  cat <<EOF
Usage: $(basename "$0") <ID_or_file> [output.pdf]

  ID_or_file   Either a svaba run ID (will find *.learn.tsv.gz files)
               or a path to a .tsv / .tsv.gz learn file directly.

  output.pdf   Optional output path. Defaults to <ID>.learn.pdf

Examples:
  $(basename "$0") my_run
  $(basename "$0") my_run.tumor.learn.tsv.gz
  $(basename "$0") my_run custom_output.pdf
EOF
  exit 1
}

[[ $# -lt 1 ]] && usage

INPUT="$1"
OUTPDF="${2:-}"

# Determine the list of learn files
FILES=()
ID=""

if [[ "$INPUT" == *.tsv.gz || "$INPUT" == *.tsv ]]; then
  # Direct file path
  [[ -f "$INPUT" ]] || { echo "ERROR: file not found: $INPUT" >&2; exit 1; }
  FILES=("$INPUT")
  # derive stem for default output name
  stem=$(basename "$INPUT")
  stem="${stem%.gz}"
  stem="${stem%.tsv}"
  stem="${stem%.learn}"
  ID="$stem"
else
  # Treat as run ID — glob for learn files
  ID="$INPUT"
  shopt -s nullglob
  FILES=( ${ID}.*.learn.tsv.gz ${ID}.learn.tsv.gz )
  shopt -u nullglob
  if [[ ${#FILES[@]} -eq 0 ]]; then
    echo "ERROR: no learn files found for ID '$ID'" >&2
    echo "  looked for: ${ID}.*.learn.tsv.gz, ${ID}.learn.tsv.gz" >&2
    exit 1
  fi
fi

[[ -z "$OUTPDF" ]] && OUTPDF="${ID}.learn.pdf"

echo "Input files: ${FILES[*]}"
echo "Output PDF:  $OUTPDF"

# Build a comma-separated list of files for R (no quoting — paths with
# commas would break, but that's exotic enough to not worry about)
FILE_LIST=$(IFS=,; echo "${FILES[*]}")

LEARN_FILES="$FILE_LIST" Rscript --vanilla - "$OUTPDF" <<'REOF'
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
outpdf <- args[1]

# The file list is passed via environment variable
file_list <- strsplit(Sys.getenv("LEARN_FILES"), ",")[[1]]

# Read and combine all files
dt <- rbindlist(lapply(file_list, function(f) {
  d <- fread(f, sep = "\t")
  # add a source column from the filename
  stem <- basename(f)
  stem <- sub("\\.learn\\.tsv\\.gz$", "", stem)
  stem <- sub("\\.learn\\.tsv$", "", stem)
  d[, source := stem]
  d
}))

cat(sprintf("Loaded %s isize observations across %d RGs from %d files\n",
            format(nrow(dt), big.mark = ","),
            uniqueN(dt$rg),
            length(file_list)))

# Summary stats per RG per source
stats <- dt[, .(
  n = .N,
  mean = mean(isize),
  median = as.double(median(isize)),
  sd = sd(isize),
  p01 = quantile(isize, 0.01),
  p99 = quantile(isize, 0.99)
), by = .(source, rg)]

cat("\nPer-RG summary:\n")
print(stats[order(source, rg)])

# Trim extreme tails for plotting (1st-99th percentile per RG)
dt_trim <- dt[, .SD[isize >= quantile(isize, 0.01) &
                     isize <= quantile(isize, 0.99)],
              by = .(source, rg)]

# Create label with n for facets
dt_trim[, facet_label := paste0(rg, " (n=", format(.N, big.mark = ","), ")"),
        by = .(source, rg)]

# Determine a reasonable number of columns for facet_wrap
n_facets <- uniqueN(dt_trim[, paste(source, rg)])
ncol_facet <- min(4, ceiling(sqrt(n_facets)))

# Plot
p <- ggplot(dt_trim, aes(x = isize, fill = source)) +
  geom_histogram(bins = 80, alpha = 0.7, color = "grey30", linewidth = 0.2) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = ncol_facet) +
  labs(
    title = "Insert size distribution per read group",
    subtitle = paste("Sources:", paste(unique(dt_trim$source), collapse = ", ")),
    x = "Insert size (bp)",
    y = "Count",
    fill = "BAM"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    strip.text = element_text(size = 7),
    legend.position = if (uniqueN(dt_trim$source) > 1) "bottom" else "none"
  )

# Size the PDF to fit the facets
n_rows <- ceiling(n_facets / ncol_facet)
pdf_h <- max(4, 1.5 + n_rows * 2.5)
pdf_w <- max(6, ncol_facet * 3)

ggsave(outpdf, p, width = pdf_w, height = pdf_h, limitsize = FALSE)
cat(sprintf("\nWrote %s (%d x %d inches)\n", outpdf, as.integer(pdf_w), as.integer(pdf_h)))
REOF
