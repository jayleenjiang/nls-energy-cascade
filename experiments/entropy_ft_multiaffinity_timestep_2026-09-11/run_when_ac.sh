#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
printf '%s\n' "$$" > "$HERE/launcher.pid"

while ! pmset -g batt | head -1 | grep -q 'AC Power'; do
  sleep 30
done

exec "$HERE/run_pipeline.sh"

