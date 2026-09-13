#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
EXPERIMENT_ROOT="${1:-$ROOT/experiments/entropy_ft_2026-08-26}"
BIN_DIR="$EXPERIMENT_ROOT/bin"
PILOT_DIR="$EXPERIMENT_ROOT/pilot"
ANALYSIS_DIR="$PILOT_DIR/analysis"
CXX="${CXX:-clang++}"
PYTHON="${PYTHON:-/opt/homebrew/bin/python3}"
EIGEN3_INCLUDE="${EIGEN3_INCLUDE:-/opt/homebrew/include/eigen3}"
LIBOMP_PREFIX="${LIBOMP_PREFIX:-/opt/homebrew/opt/libomp}"

mkdir -p "$BIN_DIR" "$PILOT_DIR" "$ANALYSIS_DIR"

"$CXX" -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I"$EIGEN3_INCLUDE" \
  -I"$LIBOMP_PREFIX/include" \
  -L"$LIBOMP_PREFIX/lib" -lomp \
  "$ROOT/flux/NLS_entropy_ft.cpp" \
  -o "$BIN_DIR/entropy_ft"

"$BIN_DIR/entropy_ft" selftest

{
  date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
  git -C "$ROOT" rev-parse HEAD
  shasum -a 256 "$ROOT/flux/NLS_entropy_ft.cpp"
  shasum -a 256 "$ROOT/flux/analyze_entropy_ft.py"
  echo 'parameters=T1=10 Tn=2 gamma=0.1 dt=0.0005 burnin=500 block_time=20 batches=8 blocks_per_stream=157 threads_per_run=2'
} > "$PILOT_DIR/pilot_manifest.txt"

pids=()
for n in 10 20 30 40; do
  prefix="$PILOT_DIR/n${n}"
  seed=$((2026082600 + n))
  "$BIN_DIR/entropy_ft" sample 10 2 "$n" 8 500 20 157 0.0005 \
    "$seed" 2 "$prefix" > "$PILOT_DIR/n${n}.log" 2>&1 &
  pids+=("$!")
done

status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
done
if [[ "$status" -ne 0 ]]; then
  echo "At least one pilot run failed; inspect $PILOT_DIR/*.log" >&2
  exit "$status"
fi

"$PYTHON" "$ROOT/flux/analyze_entropy_ft.py" \
  "$PILOT_DIR"/n*_blocks.csv \
  --output-dir "$ANALYSIS_DIR" \
  --taus 20,40,60,80,100,120,140,160,180,200 \
  --symmetric-bins 40 \
  --min-effective-count 20 \
  --bootstrap 500 \
  --seed 2026082699 \
  > "$PILOT_DIR/analysis.log" 2>&1

date -u '+completed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$PILOT_DIR/pilot_manifest.txt"
echo "Pilot complete: $PILOT_DIR"
