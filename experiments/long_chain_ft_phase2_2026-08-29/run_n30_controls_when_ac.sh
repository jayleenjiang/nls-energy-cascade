#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
runner="$repo_root/experiments/long_chain_ft_phase2_2026-08-29/run_phase2.sh"

while pmset -g batt | head -n 1 | grep -q "Battery Power"; do
  sleep 60
done

# Keep the Mac awake only after external power is available.  The simulation
# matrix and seeds remain those frozen in PROTOCOL.md.
exec caffeinate -i bash "$runner" n30-controls
