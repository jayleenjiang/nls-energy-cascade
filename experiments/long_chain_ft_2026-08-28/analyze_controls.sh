#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
root="$repo_root/experiments/long_chain_ft_2026-08-28"
aggregate="$repo_root/flux/analyze_long_chain_ft.py"
compare="$repo_root/flux/analyze_long_chain_ft_controls.py"

python3 "$aggregate" \
  "$root/selection2_n10/N1024_t80" \
  "$root/numerical_controls_n10_remedial" \
  --output-dir "$root/analysis_numerical_controls_n10" \
  --horizons 20 40 60
python3 "$compare" "$root/analysis_numerical_controls_n10" \
  --output-dir "$root/analysis_numerical_controls_n10/comparisons" \
  --horizon 60

python3 "$aggregate" \
  "$root/selection2_n40_k04_06_remedial/N512_t120" \
  "$root/numerical_controls_n40/N512_dt2p5e4_selection2" \
  "$root/numerical_controls_n40_remedial/N512_dt5e4_selection1" \
  --output-dir "$root/analysis_numerical_controls_n40" \
  --horizons 20 40 60 80 100 120
python3 "$compare" "$root/analysis_numerical_controls_n40" \
  --output-dir "$root/analysis_numerical_controls_n40/comparisons" \
  --horizon 120
