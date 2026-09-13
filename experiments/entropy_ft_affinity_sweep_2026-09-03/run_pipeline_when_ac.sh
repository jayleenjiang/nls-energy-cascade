#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-300}"

if [[ -f "$HERE/pipeline.exitcode" ]]; then
  echo "A previous pipeline exit code exists; refusing to retry automatically." >&2
  exit 76
fi

while ! pmset -g batt | head -1 | grep -q "AC Power"; do
  printf '%s waiting_for_ac\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$HERE/ac_wait.log"
  sleep "$INTERVAL_SECONDS"
done

printf '%s ac_detected\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$HERE/ac_wait.log"
exec "$HERE/run_affinity_sweep.sh"
