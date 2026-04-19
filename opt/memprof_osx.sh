#!/usr/bin/env bash
set -euo pipefail
if (( $# < 1 )); then
  echo "Usage: $0 <cmd> [args]"
  exit 1
fi
OUT_LOG="memlog.txt"
PLOT="memory.png"
: > "$OUT_LOG"
rm -f "$PLOT"
"$@" &
PID=$!
echo "[$(date)] Started PID $PID: $*" >&2
START=$(date +%s)
LAST_PLOT=$START
while kill -0 "$PID" 2>/dev/null; do
  NOW=$(date +%s)
  ELAPSED=$(( NOW - START ))
  # macOS: footprint reports physical memory in bytes (matches Activity Monitor)
MEM_MB=$(footprint -p "$PID" 2>/dev/null \
    | awk '/^[[:space:]]*phys_footprint:/ {
        val=$2; unit=$3;
        if (unit == "B")  val=val/1024/1024;
        else if (unit == "KB") val=val/1024;
        else if (unit == "GB") val=val*1024;
        # MB stays as-is
        printf "%d\n", val; exit
      }')
  MEM_MB=${MEM_MB:-0}
  # Convert to KB to match the original log format (R script divides by 1024/1024 for GB)
  MEM=$(( MEM_MB * 1024 ))  

  echo "$ELAPSED $MEM" >> "$OUT_LOG"
  if (( NOW - LAST_PLOT >= 30 )); then
    LAST_PLOT=$NOW
    if (( $(wc -l < "$OUT_LOG") >= 2 )); then
      Rscript --vanilla - <<EOF
suppressPackageStartupMessages(library(ggplot2))
df <- read.table("$OUT_LOG", col.names=c("time","mem"))
p <- ggplot(df, aes(x=time, y=mem/1024/1024)) +
     geom_line() +
     labs(x="Seconds since start", y="Physical footprint (GB)",
          title="Memory profile") +
     theme_minimal()
ggsave("$PLOT", plot=p, width=6, height=4)
EOF
      echo "[$(date)] Updated plot: $PLOT" >&2
    fi
  fi
  sleep 1
done
wait "$PID" || true
