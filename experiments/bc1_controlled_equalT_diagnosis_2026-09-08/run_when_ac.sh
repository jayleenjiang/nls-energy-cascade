#!/bin/bash
set -uo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
while ! pmset -g batt | head -1 | grep -q "AC Power"; do
  echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) waiting for AC power"
  sleep 60
done

echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) AC available; starting frozen pipeline"
"$ROOT/run_pipeline.sh"
status=$?
echo "$status" > "$ROOT/pipeline.exitcode"
echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) pipeline exit=$status"
exit "$status"

