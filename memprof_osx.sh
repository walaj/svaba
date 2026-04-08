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

  # macOS: ps -o rss gives RSS in KB
  MEM=$(ps -o rss= -p "$PID" 2>/dev/null | awk '{print $1+0}')
  MEM=${MEM:-0}

  echo "$ELAPSED $MEM" >> "$OUT_LOG"

  if (( NOW - LAST_PLOT >= 30 )); then
    LAST_PLOT=$NOW
    if (( $(wc -l < "$OUT_LOG") >= 2 )); then
      Rscript --vanilla - <<EOF
suppressPackageStartupMessages(library(ggplot2))
df <- read.table("$OUT_LOG", col.names=c("time","mem"))
p <- ggplot(df, aes(x=time, y=mem/1024/1024)) +
     geom_line() +
     labs(x="Seconds since start", y="RSS (GB)",
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
