#!/usr/bin/env bash
# svaba_explore — one-shot launcher for the bps explorer (+ optional r2c.db).
#
# Spins up a tiny localhost HTTP server in a temp directory containing
# symlinks to:
#   * docs/bps_explorer.html  (the viewer)
#   * the bps.txt[.gz] you passed
#   * the r2c.db you passed (OPTIONAL)
# …then opens the browser at
#   http://localhost:<PORT>/bps_explorer.html?bps=...[&db=...]
# The viewer's auto-loader picks the query params up and `fetch()`s
# whatever's present. With just bps you get the full filterable bps
# table, histograms, IGV links. With r2c.db too you also get the
# per-read alignment-plot side panel.
#
# Why an HTTP server (and not a self-contained file:// HTML)?
# bps.txt files are usually tiny but r2c.db can be 10s–100s of MB even
# for somatic-only filtered subsets. base64-embedding bloats by 33% AND
# the browser has to parse a multi-hundred-MB HTML file before any
# JavaScript runs. A loopback server is faster, less memory-hungry,
# and lets us stream the .db lazily.
#
# Usage:
#   svaba_explore [-p PORT] [--browser BIN] [-g GENES] <bps.txt[.gz]> [r2c.db]
#
# Examples:
#   # bps table only — works without --dump-reads outputs:
#   svaba_explore sample.bps.sorted.dedup.pass.txt.gz
#
#   # bps table + per-read alignment plot (requires --dump-reads):
#   svaba_explore sample.somatic.bps.txt.gz sample.r2c.db
#
#   svaba_explore -p 9000 sample.bps.txt sample.r2c.db
#   svaba_explore --browser 'Google Chrome' sample.bps.txt.gz
#   svaba_explore -g gencode.basic.bed.gz sample.bps.txt.gz sample.r2c.db
#
# When the wrapper is launched the URL is also printed to stdout so you
# can paste it elsewhere (e.g. into a remote-port-forwarded machine's
# Chrome). Press Ctrl+C to stop the server and clean up.

set -euo pipefail

# -------- arg parsing -------------------------------------------------------

PORT=8765
BROWSER=""    # empty → use OS default opener
BPS=""
DB=""
GENES=""      # optional .bed / .gtf gene annotation

usage() {
  cat >&2 <<'EOF'
Usage: svaba_explore [options] <bps.txt[.gz]> [r2c.db]

  <bps.txt[.gz]>        Required. The bps file to explore — typically
                        ${ID}.bps.sorted.dedup.pass.txt.gz from
                        `svaba postprocess`.
  [r2c.db]              Optional. The merged r2c SQLite database
                        produced when svaba was run with --dump-reads.
                        When provided, the viewer's per-read alignment
                        plot side-panel lights up. When omitted, you
                        still get the full filterable bps table.

  -p, --port  N         Localhost port to serve on. [8765]
      --browser BIN     Open in this browser specifically. macOS-friendly:
                        "Google Chrome", "Safari", "Firefox". Linux: pass
                        the binary name (chromium, firefox, ...). Empty
                        means OS default.
  -g, --genes FILE      Optional gene-annotation file (.bed / .bed.gz /
                        .gtf / .gtf.gz). When supplied, every variant
                        row gets a Genes column populated from this
                        file's intervals. BED is recommended (smaller +
                        faster than GTF).
  -h, --help            This message.

The bps_explorer.html source lives at <repo>/docs/bps_explorer.html
relative to this script — the wrapper resolves it automatically.
EOF
  exit 1
}

while [ $# -gt 0 ]; do
  case "$1" in
    -p|--port)    PORT="$2"; shift 2;;
    --browser)    BROWSER="$2"; shift 2;;
    -g|--genes)   GENES="$2"; shift 2;;
    -h|--help)    usage;;
    -*)           echo "Unknown option: $1" >&2; usage;;
    *)            if   [ -z "$BPS" ]; then BPS="$1"
                  elif [ -z "$DB"  ]; then DB="$1"
                  else echo "Too many positional args" >&2; usage
                  fi
                  shift;;
  esac
done

# bps is required; r2c.db is optional (when omitted, the explorer comes
# up without the per-read alignment-plot panel — bps table still works).
[ -n "$BPS" ] || usage
[ -f "$BPS" ] || { echo "Error: bps file not found: $BPS" >&2; exit 1; }
if [ -n "$DB" ]; then
  [ -f "$DB" ] || { echo "Error: r2c.db not found: $DB" >&2; exit 1; }
fi
if [ -n "$GENES" ]; then
  [ -f "$GENES" ] || { echo "Error: genes file not found: $GENES" >&2; exit 1; }
fi

# -------- locate bps_explorer.html -----------------------------------------
#
# Resolve symlinks so an installed `svaba_explore` (e.g. via `make install`)
# still finds the docs/ that ships with the repo. Falls back to a few
# typical locations if the relative path doesn't exist.

resolve_self() {
  # readlink -f isn't on stock macOS; emulate it portably.
  local src="$1"
  while [ -L "$src" ]; do
    local dir
    dir="$(cd -P "$(dirname "$src")" && pwd)"
    src="$(readlink "$src")"
    case "$src" in /*) ;; *) src="$dir/$src";; esac
  done
  echo "$(cd -P "$(dirname "$src")" && pwd)/$(basename "$src")"
}

SCRIPT="$(resolve_self "$0")"
SCRIPT_DIR="$(dirname "$SCRIPT")"

HTML_CANDIDATES=(
  "$SCRIPT_DIR/../docs/bps_explorer.html"
  "$SCRIPT_DIR/../share/svaba/bps_explorer.html"
  "$SCRIPT_DIR/bps_explorer.html"
)
HTML=""
for cand in "${HTML_CANDIDATES[@]}"; do
  if [ -f "$cand" ]; then HTML="$cand"; break; fi
done
[ -n "$HTML" ] || {
  echo "Error: cannot find bps_explorer.html. Looked in:" >&2
  printf '  %s\n' "${HTML_CANDIDATES[@]}" >&2
  exit 1
}

# Need absolute paths for symlinks to work regardless of cwd.
abspath() {
  ( cd "$(dirname "$1")" && printf '%s/%s\n' "$(pwd)" "$(basename "$1")" )
}
HTML="$(abspath "$HTML")"
BPS="$(abspath "$BPS")"
[ -n "$DB" ]    && DB="$(abspath "$DB")"
[ -n "$GENES" ] && GENES="$(abspath "$GENES")"

# -------- python3 + free-port handling -------------------------------------

command -v python3 >/dev/null 2>&1 || {
  echo "Error: python3 not in PATH (needed for the localhost http server)" >&2
  exit 1
}

# Detect a free port if the requested one is taken. Try 16 candidates
# above the requested PORT before giving up.
port_free() {
  python3 - "$1" <<'PY'
import socket, sys
port = int(sys.argv[1])
s = socket.socket()
try:
    s.bind(("127.0.0.1", port))
    print("FREE")
except OSError:
    print("BUSY")
finally:
    s.close()
PY
}

orig_port="$PORT"
for _try in $(seq 0 15); do
  if [ "$(port_free "$PORT")" = "FREE" ]; then break; fi
  PORT=$((PORT + 1))
done
if [ "$PORT" != "$orig_port" ]; then
  echo "Note: port $orig_port was busy; using $PORT" >&2
fi

# -------- temp dir + symlinks ----------------------------------------------

SERVE_DIR="$(mktemp -d -t svaba_explore.XXXXXX)"
cleanup() {
  if [ -n "${SERVER_PID:-}" ]; then
    kill "$SERVER_PID" 2>/dev/null || true
    wait "$SERVER_PID" 2>/dev/null || true
  fi
  rm -rf "$SERVE_DIR"
}
trap cleanup EXIT INT TERM

ln -s "$HTML" "$SERVE_DIR/bps_explorer.html"
ln -s "$BPS"  "$SERVE_DIR/$(basename "$BPS")"
if [ -n "$DB" ]; then
  ln -s "$DB" "$SERVE_DIR/$(basename "$DB")"
fi
if [ -n "$GENES" ]; then
  ln -s "$GENES" "$SERVE_DIR/$(basename "$GENES")"
fi

BPS_FN="$(basename "$BPS")"
url_encode() {
  python3 -c 'import sys, urllib.parse; print(urllib.parse.quote(sys.argv[1], safe=""))' "$1"
}
URL="http://localhost:${PORT}/bps_explorer.html?bps=$(url_encode "$BPS_FN")"
if [ -n "$DB" ]; then
  URL="${URL}&db=$(url_encode "$(basename "$DB")")"
fi
if [ -n "$GENES" ]; then
  URL="${URL}&genes=$(url_encode "$(basename "$GENES")")"
fi

# -------- start server + open browser --------------------------------------

cd "$SERVE_DIR"
# `--bind 127.0.0.1` keeps the server off any external interface so this
# isn't a "share my data with the LAN" situation by accident.
python3 -m http.server "$PORT" --bind 127.0.0.1 >/dev/null 2>&1 &
SERVER_PID=$!

# Tiny delay so the server has time to bind before the browser fetches.
sleep 0.3

echo "─────────────────────────────────────────────"
echo " svaba_explore"
echo "   bps:    $BPS"
if [ -n "$DB" ]; then
  echo "   db:     $DB"
else
  echo "   db:     (none — alignment-plot panel will be inactive)"
fi
[ -n "$GENES" ] && echo "   genes:  $GENES"
echo "   url:    $URL"
echo "─────────────────────────────────────────────"
echo "Browser should open automatically. If not, paste the URL above."
echo "Press Ctrl+C to stop the server."
echo

# Try to open the browser. macOS uses `open`, Linux uses xdg-open, both
# fall back to printing the URL if neither exists or --browser was given
# but isn't available.
open_url() {
  local browser="$1" url="$2" uname_s
  uname_s="$(uname -s 2>/dev/null || echo unknown)"
  if [ -n "$browser" ]; then
    case "$uname_s" in
      Darwin) open -a "$browser" "$url" 2>/dev/null && return 0 ;;
      *)      command -v "$browser" >/dev/null 2>&1 && \
              ( "$browser" "$url" >/dev/null 2>&1 & ) && return 0 ;;
    esac
  fi
  case "$uname_s" in
    Darwin) command -v open      >/dev/null 2>&1 && open "$url" && return 0 ;;
    Linux)  command -v xdg-open  >/dev/null 2>&1 && xdg-open "$url" && return 0 ;;
  esac
  return 1
}
if ! open_url "$BROWSER" "$URL"; then
  echo "(could not auto-open browser; please open the URL above manually)"
fi

# Wait for ctrl-C; trap will clean up.
wait "$SERVER_PID"
