#!/bin/bash
set -euo pipefail

ROOT=/Users/jayleenjiang/Documents/NLS
EXP="$ROOT/experiments/entropy_ft_n3_direct_2026-08-31"
SOURCE="$ROOT/flux/NLS_entropy_ft.cpp"
BIN="$EXP/bin/entropy_ft_n3"
PREFIX="$EXP/raw/n3"

mkdir -p "$EXP/bin" "$EXP/raw" "$EXP/analysis" "$EXP/provenance"

clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/include/eigen3 \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  "$SOURCE" -o "$BIN"

"$BIN" selftest | tee "$EXP/provenance/selftest.log"

COMMAND=("$BIN" sample_n3 10 2 3 8 500 20 7813 0.0005 2026083133 8 "$PREFIX" 1)
{
  date -u +started_utc=%Y-%m-%dT%H:%M:%SZ
  printf 'git_commit='; git -C "$ROOT" rev-parse HEAD
  shasum -a 256 "$SOURCE" "$BIN"
  printf 'command='; printf '%q ' "${COMMAND[@]}"; printf '\n'
  clang++ --version | head -n 1
  sysctl -n machdep.cpu.brand_string 2>/dev/null || true
} > "$EXP/provenance/production_manifest.txt"

"${COMMAND[@]}" 2>&1 | tee "$EXP/raw/n3.log"

python3 "$EXP/analyze_feasibility.py" \
  "$PREFIX"_blocks.csv "$PREFIX"_summary.csv "$EXP/analysis"

{
  date -u +completed_utc=%Y-%m-%dT%H:%M:%SZ
  shasum -a 256 "$PREFIX"_blocks.csv "$PREFIX"_summary.csv \
    "$EXP/analysis/negative_tail_counts.csv" \
    "$EXP/analysis/first_law_residuals.csv" \
    "$EXP/analysis/medium_entropy_ift.csv" \
    "$EXP/analysis/analysis_audit.json"
} >> "$EXP/provenance/production_manifest.txt"
