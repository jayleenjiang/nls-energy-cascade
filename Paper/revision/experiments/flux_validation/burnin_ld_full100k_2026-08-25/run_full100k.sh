#!/usr/bin/env bash
set -euo pipefail

cd /Users/jayleenjiang/Documents/NLS

EXP="Paper/revision/experiments/flux_validation/burnin_ld_full100k_2026-08-25"
mkdir -p "$EXP"

echo "Run started: $(date)"
echo "Experiment directory: $EXP"
echo "Host: $(hostname)"
echo "Git commit: $(git rev-parse --short HEAD 2>/dev/null || echo unknown)"
echo "Git branch: $(git branch --show-current 2>/dev/null || echo unknown)"
echo

echo "Rebuilding flux/flux_relax_tau"
clang++ -O3 -mcpu=native -std=c++17 \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/include/eigen3 \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  flux/NLS_flux_relaxation_tau.cpp -o flux/flux_relax_tau

echo
echo "Step 1: no-burn-in transient relaxation, 1024 trajectories each"
for n in 10 20 40 80; do
  prefix="$EXP/transient_n${n}"
  echo "[$(date)] transient n=$n -> $prefix"
  ./flux/flux_relax_tau transient 10 2 "$n" 64 500 0.0005 50 20260825 8 "$prefix"
done

echo
echo "Step 2: finite-tau flux-tail study, 100000 trajectories each"
for n in 10 20 30 40; do
  prefix="$EXP/tau_n${n}"
  echo "[$(date)] tau n=$n -> $prefix"
  ./flux/flux_relax_tau tau 10 2 "$n" 6250 500 200 0.0005 20 20260825 8 "$prefix"
done

echo
echo "Step 3: analyze and generate figures"
python3 flux/analyze_burnin_ld.py \
  --experiment-dir "$EXP" \
  --threshold-step 0.01 \
  --threshold-max 0.70

echo
echo "Run finished: $(date)"
