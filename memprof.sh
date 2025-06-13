#!/usr/bin/env bash
set -euo pipefail

if (( $# < 1 )); then
  echo "Usage: $0 <cmd> [args]"
  exit 1
fi

OUT_LOG="memlog.txt"
PLOT="memory.png"

# Clean slate
: > "$OUT_LOG"
rm -f "$PLOT"

# Start the users command
"$@" &
PID=$!
echo "[$(date)] Started PID $PID: $*" >&2

START=$(date +%s)
LAST_PLOT=$START

# Monitor until the process exits
while kill -0 "$PID" 2>/dev/null; do
  NOW=$(date +%s)
  ELAPSED=$(( NOW - START ))
  # Get PSS in KB
  PSS=$(awk '/^Pss:/ { print $2 }' /proc/"$PID"/smaps_rollup 2>/dev/null || echo 0)
  echo "$ELAPSED $PSS" >> "$OUT_LOG"
  
#  RSS=$(ps -o rss= -p "$PID" | awk '{print $1}')
  #echo "$ELAPSED $RSS" >> "$OUT_LOG"

  # every 60s replot
  if (( ELAPSED - (LAST_PLOT - START) >= 60 )); then
    LAST_PLOT=$NOW
    Rscript --vanilla - <<EOF
library(ggplot2)
df <- read.table("$OUT_LOG", col.names = c("time","mem"))
# plot RSS vs time
p <- ggplot(df, aes(x=time, y=mem/1024/1024)) +
     geom_line() +
     labs(x="Seconds since start", y="RSS (GB)",
          title="Memory Profile of '$*'") +
     theme_minimal()
ggsave("$PLOT", plot=p, width=6, height=4)
EOF
    echo "[$(date)] Updated plot: $PLOT" >&2
  fi

  sleep 1
done

# One final log entry (in case it just exited)
NOW=$(date +%s)
ELAPSED=$(( NOW - START ))
PSS=$(awk '/^Pss:/ { print $2 }' /proc/"$PID"/smaps_rollup 2>/dev/null || echo 0)
echo "$ELAPSED $PSS" >> "$OUT_LOG"

# Final plot
Rscript --vanilla - <<EOF
library(ggplot2)
df <- read.table("$OUT_LOG", col.names = c("time","mem"))
p <- ggplot(df, aes(x=time, y=mem/1024/1024)) +
     geom_line() +3
     labs(x="Seconds since start", y="RSS (GB)",
          title="Memory Profile of '$*'") +
     theme_minimal()
ggsave("$PLOT", plot=p, width=6, height=4)
EOF
echo "[$(date)] Final plot: $PLOT" >&2

wait "$PID"
exit $?
