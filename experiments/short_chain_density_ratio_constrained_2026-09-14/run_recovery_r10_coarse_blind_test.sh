#!/bin/bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
PY=/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11
V4_CACHE=/Users/jayleenjiang/NLS_section4_density_cache/short_chain_density_v4_holdout_2026-09-14
export TF_NUM_INTRAOP_THREADS=4 TF_NUM_INTEROP_THREADS=1 OMP_NUM_THREADS=4

if [[ -e "$ROOT/recovery_r10_coarse/driven_test" || -e "$ROOT/recovery_r10_coarse/equilibrium_test" ]]; then
  echo "refusing to reopen an existing R10 coarse blind test" >&2
  exit 2
fi

NLS_SECTION4_V4_CACHE="$V4_CACHE" "$PY" "$ROOT/scripts/evaluate_triple_calibrated_v4.py" \
  --phase driven --timestep coarse --split test \
  --model "$ROOT/recovery_r10_coarse/models/seed_8201" \
  --model "$ROOT/recovery_r10_coarse/models/seed_8202" \
  --model "$ROOT/recovery_r10_coarse/models/seed_8203" \
  --edges-from "$ROOT/recovery_r10_coarse/driven_validation/marginal_edges.npz" \
  --validation-verdict "$ROOT/recovery_r10_coarse/driven_validation/verdict.json" \
  --output "$ROOT/recovery_r10_coarse/driven_test" \
  > "$ROOT/logs/recovery_r10_coarse/driven_test.log" 2>&1

NLS_SECTION4_V4_CACHE="$V4_CACHE" "$PY" "$ROOT/scripts/evaluate_v4.py" \
  --phase equilibrium --timestep coarse --split test \
  --model "$ROOT/formal_v7_coarse/equilibrium_models/seed_10201" \
  --model "$ROOT/formal_v7_coarse/equilibrium_models/seed_10202" \
  --model "$ROOT/formal_v7_coarse/equilibrium_models/seed_10203" \
  --edges-from "$ROOT/formal_v7_coarse/equilibrium_validation/marginal_edges.npz" \
  --validation-verdict "$ROOT/formal_v7_coarse/equilibrium_validation/verdict.json" \
  --output "$ROOT/recovery_r10_coarse/equilibrium_test" \
  > "$ROOT/logs/recovery_r10_coarse/equilibrium_test.log" 2>&1
