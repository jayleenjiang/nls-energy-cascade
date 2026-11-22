#!/usr/bin/env bash
# batch_09 part B: one row of tasks_B.csv per array task (Sec. 3.3 action current, nls_cartesian_v3).
# Each task compiles its own binary on its node (-march=native, as the Sec. 3.1 runs) and runs the self test first.
# Output: results/<group>/<label>_summary.txt, <label>_current_streams.txt, <label>.log
#SBATCH --job-name=b09_B
#SBATCH --array=0-39
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=4G
#SBATCH --time=20:00:00
#SBATCH --output=logs/slurm_b09_B_%A_%a.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
CASES=tasks_B.csv
LINE=$(awk -F, -v t="$SLURM_ARRAY_TASK_ID" 'NR>1 && $1==t' "$CASES")
[ -n "$LINE" ] || { echo "no case for task $SLURM_ARRAY_TASK_ID"; exit 1; }
IFS=, read -r TASK GROUP LABEL BC T1 TN N CLOSED BURN MEASURE BATCHES DT SEED <<< "$LINE"
OUTDIR=results/$GROUP; mkdir -p "$OUTDIR"
if [ -e "$OUTDIR/${LABEL}_summary.txt" ]; then echo "refusing to overwrite $OUTDIR/${LABEL}_summary.txt"; exit 1; fi
echo "task $TASK $LABEL on $(hostname), cpus=$SLURM_CPUS_PER_TASK, $(date)"
BUILD=build/${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}; mkdir -p "$BUILD/eigen"
tar xzf eigen-3.4.0-headers.tar.gz -C "$BUILD/eigen"
g++ -O3 -march=native -fopenmp -std=c++17 -I"$BUILD/eigen/eigen3" nls_cartesian_v3.cpp -o "$BUILD/nls_cartesian"
"$BUILD/nls_cartesian" --selftest
"$BUILD/nls_cartesian" --n "$N" --T1 "$T1" --Tn "$TN" --bc "${BC#BC}" --closed "$CLOSED" --gamma 0.1 \
    --dt "$DT" --burn "$BURN" --measure "$MEASURE" --interval 0.1 --batches "$BATCHES" \
    --threads "$SLURM_CPUS_PER_TASK" --checkpoints 4 --seed "$SEED" --out "$OUTDIR/$LABEL" > "$OUTDIR/$LABEL.log" 2>&1
grep -E "midpoint_failures|nonfinite|elapsed" "$OUTDIR/${LABEL}_summary.txt"
rm -rf "$BUILD/eigen"
echo "done $LABEL $(date)"
