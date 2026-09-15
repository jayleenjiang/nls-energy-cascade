#!/bin/bash
set -euo pipefail

EXP="$(cd "$(dirname "$0")" && pwd)"
SESSION="nls_spectral_autocorr"

if tmux has-session -t "$SESSION" 2>/dev/null; then
  echo "Session $SESSION already exists" >&2
  exit 3
fi
if pgrep -f "$EXP/bin/NLS_stationary_autocorr" >/dev/null; then
  echo "A matching sampler is already running" >&2
  exit 3
fi

tmux new-session -d -s "$SESSION" \
  "/bin/bash '$EXP/run_pipeline.sh'; rc=\$?; printf '%s\n' \"\$rc\" > '$EXP/pipeline.exitcode'; exit \$rc"
tmux display-message -p -t "$SESSION" '#{session_pid}' > "$EXP/pipeline.pid"
date -u +%Y-%m-%dT%H:%M:%SZ > "$EXP/provenance/launched_utc.txt"
echo "launched $SESSION pid=$(cat "$EXP/pipeline.pid")"
