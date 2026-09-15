#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY=/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11
OLD_CACHE=/Users/jayleenjiang/NLS_section4_density_cache/short_chain_section4_2026-09-14
V4_CACHE=/Users/jayleenjiang/NLS_section4_density_cache/short_chain_density_v4_holdout_2026-09-14
V6_CACHE=/Users/jayleenjiang/NLS_section4_density_cache/short_chain_density_v6_holdout_2026-09-14
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=1 OMP_NUM_THREADS=4
mkdir -p "$ROOT/recovery_r11_fine/calibration_raw" "$ROOT/recovery_r11_fine/models" "$ROOT/logs/recovery_r11_fine"

for base in 8201 8202 8203; do
  raw="$ROOT/recovery_r11_fine/calibration_raw/seed_${base}"
  final="$ROOT/recovery_r11_fine/models/seed_${base}"
  if [[ ! -e "$raw/calibration.json" ]]; then
    NLS_SECTION4_RAW_CACHE="$OLD_CACHE" NLS_SECTION4_V4_CACHE="$V4_CACHE" \
      "$PY" "$ROOT/scripts/calibrate_three_moments_fine.py" \
      --base-model "$ROOT/recovery_r4/seed_${base}" \
      --seed "$((base+39000))" --output "$raw" \
      > "$ROOT/logs/recovery_r11_fine/calibration_seed_${base}.log" 2>&1
  fi
  if [[ ! -e "$final/model_config.json" ]]; then
    "$PY" "$ROOT/scripts/scale_triple_calibration.py" \
      --source "$raw" --factor 0.625 --output "$final"
  fi
done

if [[ ! -e "$ROOT/recovery_r11_fine/equilibrium_models/manifest.json" ]]; then
  "$PY" "$ROOT/scripts/create_zero_triple_equilibrium.py" \
    --base-model "$ROOT/formal_v4_equilibrium_fine/models/seed_9201" \
    --base-model "$ROOT/formal_v4_equilibrium_fine/models/seed_9202" \
    --base-model "$ROOT/formal_v4_equilibrium_fine/models/seed_9203" \
    --seed 49201 --seed 49202 --seed 49203 \
    --output "$ROOT/recovery_r11_fine/equilibrium_models"
fi

NLS_SECTION4_V6_CACHE="$V6_CACHE" "$PY" "$ROOT/scripts/evaluate_triple_calibrated_v6.py" \
  --phase driven --timestep fine --split validation \
  --model "$ROOT/recovery_r11_fine/models/seed_8201" \
  --model "$ROOT/recovery_r11_fine/models/seed_8202" \
  --model "$ROOT/recovery_r11_fine/models/seed_8203" \
  --output "$ROOT/recovery_r11_fine/driven_validation" \
  > "$ROOT/logs/recovery_r11_fine/driven_validation.log" 2>&1

NLS_SECTION4_V6_CACHE="$V6_CACHE" "$PY" "$ROOT/scripts/evaluate_triple_calibrated_v6.py" \
  --phase equilibrium --timestep fine --split validation \
  --model "$ROOT/recovery_r11_fine/equilibrium_models/seed_49201" \
  --model "$ROOT/recovery_r11_fine/equilibrium_models/seed_49202" \
  --model "$ROOT/recovery_r11_fine/equilibrium_models/seed_49203" \
  --output "$ROOT/recovery_r11_fine/equilibrium_validation" \
  > "$ROOT/logs/recovery_r11_fine/equilibrium_validation.log" 2>&1

