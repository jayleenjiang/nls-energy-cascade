#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPERIMENT_ROOT="${1:-$ROOT/experiments/entropy_ft_2026-08-26}"
BIN_DIR="$EXPERIMENT_ROOT/bin"
PRODUCTION_DIR="$EXPERIMENT_ROOT/production"
ANALYSIS_DIR="$PRODUCTION_DIR/analysis"
CXX="${CXX:-clang++}"
PYTHON="${PYTHON:-/opt/homebrew/bin/python3}"
EIGEN3_INCLUDE="${EIGEN3_INCLUDE:-/opt/homebrew/include/eigen3}"
LIBOMP_PREFIX="${LIBOMP_PREFIX:-/opt/homebrew/opt/libomp}"
THREADS_PER_RUN="${THREADS_PER_RUN:-2}"

BATCHES=8
BLOCKS_PER_STREAM=7813
TOTAL_STREAMS=$((BATCHES * 16))
TOTAL_BLOCKS=$((TOTAL_STREAMS * BLOCKS_PER_STREAM))

mkdir -p "$BIN_DIR" "$PRODUCTION_DIR" "$ANALYSIS_DIR"

for n in 10 20 30 40; do
  if [[ -e "$PRODUCTION_DIR/n${n}_blocks.csv" ]]; then
    echo "Refusing to overwrite existing production data: $PRODUCTION_DIR/n${n}_blocks.csv" >&2
    exit 2
  fi
done

available_kb=$(df -Pk "$EXPERIMENT_ROOT" | awk 'NR==2 {print $4}')
if [[ "$available_kb" -lt 2097152 ]]; then
  echo "At least 2 GiB of free disk space is required; found ${available_kb} KiB." >&2
  exit 3
fi

"$CXX" -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I"$EIGEN3_INCLUDE" \
  -I"$LIBOMP_PREFIX/include" \
  -L"$LIBOMP_PREFIX/lib" -lomp \
  "$ROOT/flux/NLS_entropy_ft.cpp" \
  -o "$BIN_DIR/entropy_ft_production"

"$BIN_DIR/entropy_ft_production" selftest

{
  date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
  git -C "$ROOT" rev-parse HEAD
  shasum -a 256 "$ROOT/flux/NLS_entropy_ft.cpp"
  shasum -a 256 "$ROOT/flux/analyze_entropy_ft.py"
  printf 'parameters=T1=10 Tn=2 gamma=0.1 dt=0.0005 burnin=500 block_time=20 batches=%d blocks_per_stream=%d streams=%d total_blocks_per_n=%d threads_per_run=%s\n' \
    "$BATCHES" "$BLOCKS_PER_STREAM" "$TOTAL_STREAMS" "$TOTAL_BLOCKS" "$THREADS_PER_RUN"
} > "$PRODUCTION_DIR/production_manifest.txt"

pids=()
for n in 10 20 30 40; do
  prefix="$PRODUCTION_DIR/n${n}"
  seed=$((2026082700 + n))
  "$BIN_DIR/entropy_ft_production" sample 10 2 "$n" "$BATCHES" 500 20 \
    "$BLOCKS_PER_STREAM" 0.0005 "$seed" "$THREADS_PER_RUN" "$prefix" \
    > "$PRODUCTION_DIR/n${n}.log" 2>&1 &
  pids+=("$!")
  printf 'n=%d pid=%d\n' "$n" "$!" >> "$PRODUCTION_DIR/production_pids.txt"
done

status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
done
if [[ "$status" -ne 0 ]]; then
  date -u '+failed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$PRODUCTION_DIR/production_manifest.txt"
  echo "At least one production run failed; inspect $PRODUCTION_DIR/*.log" >&2
  exit "$status"
fi

"$PYTHON" "$ROOT/flux/analyze_entropy_ft.py" \
  "$PRODUCTION_DIR"/n*_blocks.csv \
  --output-dir "$ANALYSIS_DIR" \
  --taus 20,40,60,80,100,120,140,160,180,200 \
  --symmetric-bins 60 \
  --min-effective-count 50 \
  --bootstrap 1000 \
  --seed 2026082799 \
  --threshold-step 0.01 \
  > "$PRODUCTION_DIR/analysis.log" 2>&1

date -u '+completed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$PRODUCTION_DIR/production_manifest.txt"
echo "Production complete: $PRODUCTION_DIR"
