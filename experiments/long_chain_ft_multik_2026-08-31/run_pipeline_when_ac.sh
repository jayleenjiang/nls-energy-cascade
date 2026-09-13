#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
phase3="$repo_root/experiments/long_chain_ft_multik_2026-08-31"
lock="$phase3/.pipeline_lock"

if ! mkdir "$lock" 2>/dev/null; then
  echo "pipeline lock exists; refusing duplicate launch" >&2
  exit 1
fi
trap 'rmdir "$lock" 2>/dev/null || true' EXIT

while pmset -g batt | head -n 1 | grep -q "Battery Power"; do
  sleep 60
done

date -u +'%Y-%m-%dT%H:%M:%SZ' > "$phase3/STARTED_UTC.txt"

bash "$repo_root/experiments/long_chain_ft_phase2_2026-08-29/build_and_selftest.sh"
bash "$phase3/run_multik.sh" pilot
python3 "$phase3/freeze_multik_grid.py"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$phase3/PRODUCTION_STARTED"
bash "$phase3/run_multik.sh" production
bash "$phase3/run_multik.sh" population
bash "$phase3/run_multik.sh" n40-controls
bash "$phase3/analyze_multik.sh"
date -u +'%Y-%m-%dT%H:%M:%SZ' > "$phase3/COMPLETED_UTC.txt"

