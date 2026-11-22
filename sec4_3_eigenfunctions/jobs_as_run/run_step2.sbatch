#!/usr/bin/env bash
# Section 4.4 Step 2: streaming gEDMD production (n=3), one SLURM array task per case.
# Submit from this folder:  sbatch run_step2.sbatch
#SBATCH --job-name=n3gedmd
#SBATCH --array=0-3
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=slurm_%A_%a.log
set -euo pipefail
NCPU=${SLURM_CPUS_PER_TASK:-16}
EIGEN_INC=${EIGEN_INC:-/usr/include/eigen3}
COMMON="--burn 50 --measure 10000 --snapshot 0.1 --batches 16 --groups 4 --d 5 --M 3 --gamma 0.1 --threads $NCPU"
CASES=(
  "neq_dt1e-3|--T1 2 --T3 8 --dt 1e-3   --seed 5101 --s1 1.32 --s2 1.71 --s3 1.65"
  "neq_dt2p5e-4|--T1 2 --T3 8 --dt 2.5e-4 --seed 5102 --s1 1.32 --s2 1.71 --s3 1.65"
  "eq_dt1e-3|--T1 5 --T3 5 --dt 1e-3   --seed 5103 --s1 1.48 --s2 1.68 --s3 1.48"
  "eq_dt2p5e-4|--T1 5 --T3 5 --dt 2.5e-4 --seed 5104 --s1 1.48 --s2 1.68 --s3 1.48"
)
C=${CASES[$SLURM_ARRAY_TASK_ID]}
OUT=${C%%|*}; ARGS=${C#*|}
BIN=gedmd_stream_${SLURM_ARRAY_TASK_ID}
echo "task $SLURM_ARRAY_TASK_ID -> $OUT on $(hostname), cpus=$NCPU, $(date)"
g++ -O3 -march=native -fopenmp -std=c++17 -I"$EIGEN_INC" NLS3_gedmd_stream.cpp -o "$BIN"
# shellcheck disable=SC2086
./"$BIN" --out "$OUT" $ARGS $COMMON > "$OUT.log" 2>&1
tail -1 "$OUT.log"
python3 analyze_step2.py "$OUT" > "${OUT}_analysis.log" 2>&1
cat "${OUT}_analysis.log"
echo "done $OUT $(date)"
