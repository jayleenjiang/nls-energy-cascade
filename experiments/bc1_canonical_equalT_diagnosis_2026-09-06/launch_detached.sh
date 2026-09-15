#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
if [[ -f "$ROOT/launcher.pid" ]] && kill -0 "$(cat "$ROOT/launcher.pid")" 2>/dev/null; then
  echo "launcher already active: $(cat "$ROOT/launcher.pid")"
  exit 0
fi
nohup "$ROOT/run_when_ac.sh" > "$ROOT/launcher.log" 2>&1 &
echo $! > "$ROOT/launcher.pid"
echo "launched PID $(cat "$ROOT/launcher.pid")"

