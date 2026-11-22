#!/usr/bin/env bash
# Long independent validation trajectories (nls_cartesian v2): 256 streams, burn-in 500, measurement 10000
# sampled every 0.5; threads = 1 so the dump order is (checkpoint, batch, snapshot, lane).
#SBATCH --job-name=v5val
#SBATCH --array=0-7
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=03:00:00
#SBATCH --output=logs/slurm_val_%A_%a.log
set -euo pipefail
CASES=("3 5 5 eq A 1" "3 5 5 eq B 2" "3 2 8 neq A 3" "3 2 8 neq B 4" "4 5 5 eq A 1" "4 5 5 eq B 2" "4 2 8 neq A 3" "4 2 8 neq B 4")
read -r N T1 TN CASE SET K <<< "${CASES[$SLURM_ARRAY_TASK_ID]}"
OUT=val_long/n${N}_${CASE}_val${SET}
mkdir -p val_long
if [ -e "${OUT}_summary.txt" ]; then echo "refusing to overwrite $OUT"; exit 1; fi
BUILD=build/val_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p "$BUILD/eigen"
tar xzf eigen-3.4.0-headers.tar.gz -C "$BUILD/eigen"
g++ -O3 -march=native -fopenmp -std=c++17 -I"$BUILD/eigen/eigen3" nls_cartesian.cpp -o "$BUILD/nls_cartesian"
"$BUILD/nls_cartesian" --selftest
"$BUILD/nls_cartesian" --n "$N" --T1 "$T1" --Tn "$TN" --bc 1 --dt 1e-3 --burn 500 --measure 10000 --interval 0.5 \
    --batches 16 --threads 1 --checkpoints 4 --dump 1 --seed 2026101150${N}${K} --out "$OUT" > "$OUT.log" 2>&1
grep -E "midpoint_failures|nonfinite|elapsed|^mean_I" "${OUT}_summary.txt"
rm -rf "$BUILD/eigen"
