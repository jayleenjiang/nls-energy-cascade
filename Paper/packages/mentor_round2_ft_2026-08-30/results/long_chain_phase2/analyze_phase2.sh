#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
phase1="$repo_root/experiments/long_chain_ft_2026-08-28"
phase2="$repo_root/experiments/long_chain_ft_phase2_2026-08-29"
raw="$phase2/raw"
analysis="$phase2/analysis/final"

analyze_n20() {
  python3 "$repo_root/flux/analyze_long_chain_ft.py" \
    "$phase1/selection2_n20_k04_06/N1024_t60" \
    "$raw/n20_dt_N1024" "$raw/n20_selection1_N1024" \
    --output-dir "$analysis/n20_controls" --horizons 20 40 60
  python3 "$repo_root/flux/analyze_long_chain_ft_controls.py" \
    "$analysis/n20_controls" --output-dir "$analysis/n20_controls/comparisons" \
    --horizon 60 --baseline-dt 0.0005 --baseline-selection 2
}

analyze_n30_population() {
  python3 "$repo_root/flux/analyze_long_chain_ft.py" \
    "$phase1/selection2_n30_k04_06/N1024_t60" \
    "$phase1/selection2_n30_k04_06_remedial/N2048_t60" \
    "$raw/n30_population_N4096" \
    --output-dir "$analysis/n30_population" --horizons 20 40 60
}

analyze_n30_controls() {
  python3 "$repo_root/flux/analyze_long_chain_ft.py" \
    "$raw/n30_population_N4096" "$raw/n30_dt_N4096" \
    "$raw/n30_selection1_N4096" \
    --output-dir "$analysis/n30_controls" --horizons 20 40 60
  python3 "$repo_root/flux/analyze_long_chain_ft_controls.py" \
    "$analysis/n30_controls" --output-dir "$analysis/n30_controls/comparisons" \
    --horizon 60 --baseline-dt 0.0005 --baseline-selection 2
}

analyze_n40() {
  python3 "$repo_root/flux/analyze_long_chain_ft.py" \
    "$phase1/selection2_n40_k04_06_remedial/N1024_t120" \
    "$raw/n40_selection1_N1024" "$raw/n40_dt_N1024" \
    --output-dir "$analysis/n40_controls" --horizons 20 40 60 80 100 120
  python3 "$repo_root/flux/analyze_long_chain_ft_controls.py" \
    "$analysis/n40_controls" --output-dir "$analysis/n40_controls/comparisons" \
    --horizon 120 --baseline-dt 0.0005 --baseline-selection 2
}

target=${1:-}
case "$target" in
  stage1)
    analyze_n20
    analyze_n30_population
    analyze_n40
    ;;
  final)
    analyze_n20
    analyze_n30_population
    analyze_n30_controls
    analyze_n40
    ;;
  *)
    echo "usage: $0 {stage1|final}" >&2
    exit 2
    ;;
esac
