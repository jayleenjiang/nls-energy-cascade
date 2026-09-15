#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
root="$repo_root/experiments/long_chain_ft_2026-08-28"
analyzer="$repo_root/flux/analyze_long_chain_ft.py"

python3 "$analyzer" "$root/selection2_n10" \
  --output-dir "$root/analysis_selection2_n10_population" \
  --horizons 20 40 60

python3 "$analyzer" "$root/selection2_n20_k04_06" \
  --output-dir "$root/analysis_n20_k04_06_population" \
  --horizons 20 40 60

python3 "$analyzer" "$root/selection2_n30_k04_06" \
  --output-dir "$root/analysis_n30_k04_06_population" \
  --horizons 20 40 60

if [[ -d "$root/selection2_n30_k04_06_remedial/N2048_t60" ]]; then
  python3 "$analyzer" \
    "$root/selection2_n30_k04_06/N1024_t60" \
    "$root/selection2_n30_k04_06_remedial/N2048_t60" \
    --output-dir "$root/analysis_n30_k04_06_remedial" \
    --horizons 20 40 60
fi

python3 "$analyzer" \
  "$root/selection2_n40_k04_06/N512_t60" \
  "$root/selection2_n40_k04_06/N1024_t80" \
  --output-dir "$root/analysis_n40_k04_06_population_t60" \
  --horizons 20 40 60

# The extended audit keeps the t=80 N_c=1024 result visible without replacing
# the common t=60 population comparison.
python3 "$analyzer" \
  "$root/selection2_n40_k04_06/N512_t60" \
  "$root/selection2_n40_k04_06/N1024_t80" \
  --output-dir "$root/analysis_n40_k04_06_extended" \
  --horizons 20 40 60 80

if [[ -d "$root/selection2_n40_k04_06_remedial/N512_t120" && \
      -d "$root/selection2_n40_k04_06_remedial/N1024_t120" ]]; then
  python3 "$analyzer" "$root/selection2_n40_k04_06_remedial" \
    --output-dir "$root/analysis_n40_k04_06_remedial" \
    --horizons 20 40 60 80 100 120
fi
