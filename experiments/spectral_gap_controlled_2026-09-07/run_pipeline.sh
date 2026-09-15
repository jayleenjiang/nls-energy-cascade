#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
SRC="$ROOT/source/NLS_stationary_autocorr_controlled.cpp"
BIN="$ROOT/bin/NLS_stationary_autocorr_controlled"
mkdir -p "$ROOT/bin" "$ROOT/raw" "$ROOT/analysis" "$ROOT/provenance"

clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/include/eigen3 \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  "$SRC" -o "$BIN"

shasum -a 256 "$SRC" "$BIN" > "$ROOT/provenance/source_binary.sha256"
git -C "$ROOT" rev-parse HEAD > "$ROOT/provenance/launch_commit.txt"
date -u +%Y-%m-%dT%H:%M:%SZ > "$ROOT/provenance/launched_utc.txt"

run_case() {
  local name="$1" T1="$2" T3="$3" dt="$4" seed="$5"
  local out="$ROOT/raw/$name"
  if [[ -e "$out/metadata.csv" || -e "$out/stream_000.f32" ]]; then
    echo "refusing to overwrite existing production output: $out" >&2
    return 3
  fi
  mkdir -p "$out"
  echo "$BIN $out $T1 $T3 0.1 $dt 500 1000 0.01 4 $seed 4" \
    | tee -a "$ROOT/provenance/production_commands.txt"
  "$BIN" "$out" "$T1" "$T3" 0.1 "$dt" 500 1000 0.01 4 "$seed" 4 \
    2>&1 | tee "$out/run.log"
}

run_case driven_dt1e-3 2 8 0.001 2026090801
run_case driven_dt2p5e-4 2 8 0.00025 2026090802
run_case equilibrium_dt1e-3 5 5 0.001 2026090803
run_case equilibrium_dt2p5e-4 5 5 0.00025 2026090804

python3 "$ROOT/analyze_controlled.py" 2>&1 | tee "$ROOT/analysis/analysis.log"
date -u +%Y-%m-%dT%H:%M:%SZ > "$ROOT/provenance/completed_utc.txt"
