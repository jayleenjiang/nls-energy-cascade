#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 PID" >&2
  exit 2
fi

PID="$1"
while kill -0 "$PID" 2>/dev/null; do
  if pmset -g batt | head -n 1 | grep -q "AC Power"; then
    kill -CONT "$PID"
    echo "resumed PID $PID on AC at $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    exit 0
  fi
  sleep 60
done

echo "PID $PID exited before AC became available" >&2
exit 1
