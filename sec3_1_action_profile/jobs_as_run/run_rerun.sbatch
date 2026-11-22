#!/usr/bin/env bash
# Section 3.1: rerun of cases that failed with the original binary, using nls_cartesian v2 + the midpoint
# non-convergence safeguard (step halving down to dt/64, exit on non-finite state).  Same seeds and settings.
#   atlas   : open chain, BC1/BC2/BC3/BC3b x (T1,Tn) in {(10,5),(4,6),(6,6)} x n in {25,50,100},
#             two seeded replicates of 16 batches (256 streams) each -> 512 trajectories
#   closed  : one-bath closed chain (left BC1 bath, right end free), T_left in {6,10}, same n / replicates
#   dense   : (4,6), BC1-BC3b, n = 20,30,...,150, burn-in 3.2 n^2, 64 trajectories
#   dtcheck : (10,5), n = 100, dt = 5e-4, 64 trajectories (time-step check against the dt = 1e-3 atlas)
# gamma = 0.1, dt = 1e-3 (except dtcheck), measurement 2000, sampling interval 0.1, 4 checkpoints.
# Each task compiles its own binary on its node (-march=native) and runs the self test first.
#SBATCH --job-name=s31_rerun3
#SBATCH --array=0-1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=4G
#SBATCH --time=08:00:00
#SBATCH --output=logs/slurm_s31_rerun3_%A_%a.log
set -euo pipefail
CASES=cases_rerun.csv
LINE=$(awk -F, -v t="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==t' "$CASES")
[ -n "$LINE" ] || { echo "no case for task $SLURM_ARRAY_TASK_ID"; exit 1; }
IFS=, read -r TASK GROUP LABEL BC T1 TN N CLOSED BURN MEASURE BATCHES DT SEED RUNM <<< "$LINE"
OUTDIR=results/$GROUP; mkdir -p "$OUTDIR"
if [ -e "$OUTDIR/${LABEL}_summary.txt" ]; then echo "refusing to overwrite $OUTDIR/${LABEL}_summary.txt"; exit 1; fi
echo "task $TASK $LABEL on $(hostname), cpus=$SLURM_CPUS_PER_TASK, $(date)"
BUILD=build/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p "$BUILD/eigen"
tar xzf eigen-3.4.0-headers.tar.gz -C "$BUILD/eigen"
g++ -O3 -march=native -fopenmp -std=c++17 -I"$BUILD/eigen/eigen3" nls_cartesian.cpp -o "$BUILD/nls_cartesian"
"$BUILD/nls_cartesian" --selftest
"$BUILD/nls_cartesian" --n "$N" --T1 "$T1" --Tn "$TN" --bc "${BC#BC}" --closed "$CLOSED" --gamma 0.1 \
    --dt "$DT" --burn "$BURN" --measure "$MEASURE" --interval 0.1 --batches "$BATCHES" \
    --threads "$SLURM_CPUS_PER_TASK" --checkpoints 4 --seed "$SEED" --runaway-M "$RUNM" --out "$OUTDIR/$LABEL" > "$OUTDIR/$LABEL.log" 2>&1
grep -E "midpoint_failures|midpoint_split_steps|nonfinite|elapsed|runaway|streams_used" "$OUTDIR/${LABEL}_summary.txt"
rm -rf "$BUILD/eigen"
echo "done $LABEL $(date)"
