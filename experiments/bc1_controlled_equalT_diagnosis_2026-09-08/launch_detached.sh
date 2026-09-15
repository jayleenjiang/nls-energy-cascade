#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
if [[ -s "$ROOT/launcher.pid" ]]; then
  old_pid="$(cat "$ROOT/launcher.pid")"
  if kill -0 "$old_pid" 2>/dev/null; then
    echo "launcher already active: $old_pid" >&2
    exit 3
  fi
fi
if pgrep -f "$ROOT/run_pipeline.sh|NLS_bc1_cartesian_profile.*$ROOT/raw" >/dev/null; then
  echo "matching diagnostic process already active" >&2
  exit 4
fi

nohup caffeinate -dimsu "$ROOT/run_when_ac.sh" \
  >> "$ROOT/launcher.log" 2>&1 < /dev/null &
echo "$!" > "$ROOT/launcher.pid"
echo "launched PID $!"
