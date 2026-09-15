#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
phase3="$repo_root/experiments/long_chain_ft_multik_2026-08-31"
raw="$phase3/raw"
analysis="$phase3/analysis"

[[ -s "$phase3/FROZEN_GRID.csv" ]] || {
  echo "FROZEN_GRID.csv is required." >&2
  exit 1
}

mkdir -p "$analysis"
python3 "$repo_root/flux/analyze_long_chain_ft.py" \
  "$repo_root/experiments/long_chain_ft_2026-08-28/selection2_n10/N2048_t60" \
  "$repo_root/experiments/long_chain_ft_2026-08-28/selection2_n20_k04_06/N1024_t60" \
  "$repo_root/experiments/long_chain_ft_phase2_2026-08-29/raw/n30_population_N4096" \
  "$repo_root/experiments/long_chain_ft_2026-08-28/selection2_n40_k04_06_remedial/N1024_t120" \
  "$raw/production" "$raw/population" \
  --output-dir "$analysis/primary" --horizons 20 40 60 80 100 120

n40_outer=$(find "$raw/production" -mindepth 1 -maxdepth 1 -type d \
  -name 'n40_outer_*' -print -quit)
[[ -n $n40_outer ]] || {
  echo "missing n=40 outer production directory" >&2
  exit 1
}
python3 "$repo_root/flux/analyze_long_chain_ft.py" \
  "$n40_outer" "$raw/n40_controls" \
  --output-dir "$analysis/n40_outer_controls" \
  --horizons 20 40 60 80 100 120
python3 "$repo_root/flux/analyze_long_chain_ft_controls.py" \
  "$analysis/n40_outer_controls" \
  --output-dir "$analysis/n40_outer_comparisons" \
  --horizon 120 --baseline-dt 0.0005 --baseline-selection 2

python3 "$phase3/validate_multik.py"

