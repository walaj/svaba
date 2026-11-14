#!/usr/bin/env Rscript
# Plot svaba thread-scaling results produced by opt/thread_scaling.sh.
#
# Usage: Rscript opt/thread_scaling_plot.R <thread_scaling.tsv> <outdir> [physical_cores_per_socket] [sockets]
#
# Writes to <outdir>:
#   thread_scaling.png   4-panel figure (wall, CPU, speedup, efficiency)
#   thread_scaling_summary.tsv  table with derived metrics
suppressPackageStartupMessages({ library(ggplot2); library(data.table) })

args <- commandArgs(trailingOnly = TRUE)
tsv    <- if (length(args) >= 1) args[1] else "/mnt/wgs/svaba_bench/thread_scaling.tsv"
outdir <- if (length(args) >= 2) args[2] else "opt/thread_scaling"
cps    <- if (length(args) >= 3) as.numeric(args[3]) else 8   # physical cores / socket
nsock  <- if (length(args) >= 4) as.numeric(args[4]) else 2
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

d <- fread(tsv)[exit == 0][order(threads)]
d <- d[, .SD[.N], by = threads]                 # last successful run per thread count
base <- d[1]                                    # smallest thread count = baseline
d[, `:=`(
  wall_min    = wall_s / 60,
  cpu_hr      = cpu_s / 3600,
  speedup     = base$wall_s / wall_s,           # relative to baseline run
  ideal       = threads / base$threads,
  efficiency  = (base$wall_s * base$threads) / (wall_s * threads),
  cpu_per_thr = cpu_s / base$cpu_s,             # CPU-time inflation vs baseline
  utilization = pct_cpu / (100 * threads),      # busy fraction of requested threads
  max_rss_gb  = max_rss_kb / 1024^2
)]
fwrite(d, file.path(outdir, "thread_scaling_summary.tsv"), sep = "\t")

# Topology markers: one socket's physical cores, all physical cores.
vlines <- data.frame(x = c(cps, cps * nsock),
                     lab = c(sprintf("%d (1 socket)", cps), sprintf("%d (all phys. cores)", cps * nsock)))
topo <- list(
  geom_vline(data = vlines, aes(xintercept = x), linetype = "dashed", colour = "grey55"),
  geom_text(data = vlines, aes(x = x, y = Inf, label = lab), angle = 90, vjust = -0.4,
            hjust = 1.1, size = 3, colour = "grey40"),
  scale_x_continuous(breaks = d$threads),
  theme_bw(base_size = 11), theme(panel.grid.minor = element_blank())
)
pt <- list(geom_line(colour = "#2b6cb0", linewidth = 0.9), geom_point(colour = "#2b6cb0", size = 2))

p1 <- ggplot(d, aes(threads, wall_min)) + pt + topo +
  labs(title = "Wall time", x = "threads (-p)", y = "minutes")
p2 <- ggplot(d, aes(threads, cpu_hr)) + pt + topo +
  labs(title = "Total CPU time (user + sys)", x = "threads (-p)", y = "CPU-hours")
p3 <- ggplot(d, aes(threads)) +
  geom_line(aes(y = ideal), linetype = "dotted", colour = "grey30") +
  geom_line(aes(y = speedup), colour = "#2b6cb0", linewidth = 0.9) +
  geom_point(aes(y = speedup), colour = "#2b6cb0", size = 2) + topo +
  labs(title = sprintf("Speedup vs -p %d (dotted = ideal)", base$threads),
       x = "threads (-p)", y = "speedup")
p4 <- ggplot(d, aes(threads, efficiency)) + pt + topo +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  labs(title = "Parallel efficiency", x = "threads (-p)",
       y = sprintf("efficiency vs -p %d", base$threads))

# Assemble 2x2 without extra packages.
png(file.path(outdir, "thread_scaling.png"), width = 12, height = 9, units = "in", res = 150)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(2, 2)))
vp <- function(r, c) grid::viewport(layout.pos.row = r, layout.pos.col = c)
print(p1, vp = vp(1, 1)); print(p2, vp = vp(1, 2))
print(p3, vp = vp(2, 1)); print(p4, vp = vp(2, 2))
invisible(dev.off())

p5 <- ggplot(d, aes(threads, max_rss_gb)) + pt + topo +
  labs(title = "Peak memory (max RSS)", x = "threads (-p)", y = "GiB")
ggsave(file.path(outdir, "thread_scaling_memory.png"), p5, width = 7, height = 4.5, dpi = 150)

print(d[, .(threads, wall_min = round(wall_min, 1), cpu_hr = round(cpu_hr, 2),
            speedup = round(speedup, 2), efficiency = round(efficiency, 2),
            utilization = round(utilization, 2), max_rss_gb = round(max_rss_gb, 1), n_bps)])
