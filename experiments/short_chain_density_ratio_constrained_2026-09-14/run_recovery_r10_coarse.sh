#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY=/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11
CACHE=/Users/jayleenjiang/NLS_section4_density_cache/short_chain_section4_2026-09-14
V4_CACHE=/Users/jayleenjiang/NLS_section4_density_cache/short_chain_density_v4_holdout_2026-09-14
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=1 OMP_NUM_THREADS=4
mkdir -p "$ROOT/recovery_r10_coarse/calibration_raw" "$ROOT/recovery_r10_coarse/models" "$ROOT/logs/recovery_r10_coarse"

for base in 8201 8202 8203; do
  raw="$ROOT/recovery_r10_coarse/calibration_raw/seed_${base}"
  final="$ROOT/recovery_r10_coarse/models/seed_${base}"
  if [[ ! -e "$raw/calibration.json" ]]; then
    NLS_SECTION4_RAW_CACHE="$CACHE" "$PY" "$ROOT/scripts/calibrate_three_moments_coarse.py" \
      --base-model "$ROOT/formal_v7_coarse/r4_base/seed_${base}" \
      --seed "$((base+38000))" --output "$raw" \
      > "$ROOT/logs/recovery_r10_coarse/calibration_seed_${base}.log" 2>&1
  fi
  if [[ ! -e "$final/model_config.json" ]]; then
    "$PY" "$ROOT/scripts/scale_triple_calibration.py" \
      --source "$raw" --factor 0.625 --output "$final"
  fi
done

NLS_SECTION4_V4_CACHE="$V4_CACHE" "$PY" "$ROOT/scripts/evaluate_triple_calibrated_v4.py" \
  --phase driven --timestep coarse --split validation \
  --model "$ROOT/recovery_r10_coarse/models/seed_8201" \
  --model "$ROOT/recovery_r10_coarse/models/seed_8202" \
  --model "$ROOT/recovery_r10_coarse/models/seed_8203" \
  --output "$ROOT/recovery_r10_coarse/driven_validation" \
  > "$ROOT/logs/recovery_r10_coarse/driven_validation.log" 2>&1

cp -R "$ROOT/formal_v7_coarse/equilibrium_validation" "$ROOT/recovery_r10_coarse/equilibrium_validation"
