#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
INTERVAL_SECONDS="${INTERVAL_SECONDS:-300}"

if [[ -f "$HERE/equilibrium_pipeline.exitcode" ]]; then
  echo "An equilibrium pipeline exit code already exists; refusing to retry." >&2
  exit 76
fi

if ! mkdir "$HERE/.equilibrium.queue.lock" 2>/dev/null; then
  echo "An equilibrium queue already exists." >&2
  exit 77
fi

printf '%s\n' "$$" > "$HERE/equilibrium_queue.pid"

while [[ ! -f "$HERE/pipeline.exitcode" ]]; do
  printf '%s waiting_for_driven_pipeline\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$HERE/equilibrium_queue.log"
  sleep "$INTERVAL_SECONDS"
done

driven_status="$(tr -d '[:space:]' < "$HERE/pipeline.exitcode")"
if [[ "$driven_status" != "0" ]]; then
  echo "Driven pipeline exited with status $driven_status; equilibrium remains queued but was not started." >&2
  exit 78
fi

while ! pmset -g batt | head -1 | grep -q "AC Power"; do
  printf '%s waiting_for_ac\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$HERE/equilibrium_queue.log"
  sleep "$INTERVAL_SECONDS"
done

printf '%s starting_equilibrium\n' "$(date -u '+%Y-%m-%dT%H:%M:%SZ')" >> "$HERE/equilibrium_queue.log"
exec "$HERE/run_equilibrium_n10.sh"
