#!/usr/bin/env bash
# Section 3.2 (local equilibrium) rerun with nls_cartesian v2, one row of cases_lte.csv per array task.
# Protocol of fixed_reference_2026-09-18: BC1 open chain, NESS (10,5) vs equilibrium (7.5,7.5),
# n = 25, 50, 100, left sites floor(0.24 n), floor(0.48 n), floor(0.72 n), burn-in 2000/8000/32000,
# measurement 2000 sampled every 0.01, two seeded replicates x 16 batches x 16 streams, one 80^3
# histogram per batch and site.  Integrator: v2 (co-rotating implicit midpoint + Euler-Maruyama baths),
# dt = 1e-3, gamma = 0.1.
#SBATCH --job-name=lte_v2
#SBATCH --array=0-11
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=logs/slurm_lte_%A_%a.log
set -euo pipefail
LINE=$(awk -F, -v t="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==t' cases_lte.csv)
[ -n "$LINE" ] || { echo "no case for task $SLURM_ARRAY_TASK_ID"; exit 1; }
IFS=, read -r TASK TEMP T1 TN N REP BATCHES BURN MEASURE INTERVAL SEED SITES <<< "$LINE"
OUTDIR=raw/cartesian_v2; mkdir -p "$OUTDIR"
PREFIX=$OUTDIR/${TEMP}_n${N}_rep${REP}
if [ -e "${PREFIX}_summary.txt" ]; then echo "refusing to overwrite ${PREFIX}_summary.txt"; exit 1; fi
echo "task $TASK $PREFIX on $(hostname), cpus=$SLURM_CPUS_PER_TASK, $(date)"
BUILD=build/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p "$BUILD/eigen"
tar xzf eigen-3.4.0-headers.tar.gz -C "$BUILD/eigen"
g++ -O3 -march=native -fopenmp -std=c++17 -I"$BUILD/eigen/eigen3" nls_cartesian.cpp -o "$BUILD/nls_cartesian"
"$BUILD/nls_cartesian" --selftest
"$BUILD/nls_cartesian" --n "$N" --T1 "$T1" --Tn "$TN" --bc 1 --gamma 0.1 --dt 1e-3 \
    --burn "$BURN" --measure "$MEASURE" --interval "$INTERVAL" --batches "$BATCHES" \
    --threads "$SLURM_CPUS_PER_TASK" --checkpoints 4 --seed "$SEED" \
    --lte-sites "${SITES// /,}" --out "$PREFIX" > "$PREFIX.log" 2>&1
grep -E "midpoint_failures|nonfinite|elapsed|lte_" "${PREFIX}_summary.txt"
rm -rf "$BUILD/eigen"
echo "done $PREFIX $(date)"
