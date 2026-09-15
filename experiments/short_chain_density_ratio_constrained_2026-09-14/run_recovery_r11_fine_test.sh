#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY=/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11
V6_CACHE=/Users/jayleenjiang/NLS_section4_density_cache/short_chain_density_v6_holdout_2026-09-14
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=1 OMP_NUM_THREADS=4
if [[ -e "$ROOT/recovery_r11_fine/driven_test" || -e "$ROOT/recovery_r11_fine/equilibrium_test" ]]; then
  echo "refusing to reopen the declared final R11 test" >&2
  exit 2
fi

NLS_SECTION4_V6_CACHE="$V6_CACHE" "$PY" "$ROOT/scripts/evaluate_triple_calibrated_v6.py" \
  --phase driven --timestep fine --split test \
  --model "$ROOT/recovery_r11_fine/models/seed_8201" \
  --model "$ROOT/recovery_r11_fine/models/seed_8202" \
  --model "$ROOT/recovery_r11_fine/models/seed_8203" \
  --edges-from "$ROOT/recovery_r11_fine/driven_validation/marginal_edges.npz" \
  --validation-verdict "$ROOT/recovery_r11_fine/driven_validation/verdict.json" \
  --output "$ROOT/recovery_r11_fine/driven_test" \
  > "$ROOT/logs/recovery_r11_fine/driven_test.log" 2>&1

NLS_SECTION4_V6_CACHE="$V6_CACHE" "$PY" "$ROOT/scripts/evaluate_triple_calibrated_v6.py" \
  --phase equilibrium --timestep fine --split test \
  --model "$ROOT/recovery_r11_fine/equilibrium_models/seed_49201" \
  --model "$ROOT/recovery_r11_fine/equilibrium_models/seed_49202" \
  --model "$ROOT/recovery_r11_fine/equilibrium_models/seed_49203" \
  --edges-from "$ROOT/recovery_r11_fine/equilibrium_validation/marginal_edges.npz" \
  --validation-verdict "$ROOT/recovery_r11_fine/equilibrium_validation/verdict.json" \
  --output "$ROOT/recovery_r11_fine/equilibrium_test" \
  > "$ROOT/logs/recovery_r11_fine/equilibrium_test.log" 2>&1

