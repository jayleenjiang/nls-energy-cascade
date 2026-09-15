#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
printf '%s\n' "$$" > "$HERE/recovery_launcher.pid"

while ! pmset -g batt | head -1 | grep -q 'AC Power'; do
  sleep 30
done

exec "$HERE/run_missing_fine_pipeline.sh"
