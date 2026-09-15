#!/bin/bash
set -euo pipefail

EXP="$(cd "$(dirname "$0")" && pwd)"
SRC="$EXP/source/NLS_stationary_autocorr.cpp"
BIN="$EXP/bin/NLS_stationary_autocorr"

if find "$EXP/raw/driven_T2_T8" "$EXP/raw/equilibrium_T5_T5" \
    -name 'stream_*.f32' -print -quit 2>/dev/null | grep -q .; then
  echo "Refusing to overwrite existing production streams" >&2
  exit 3
fi

clang++ -O3 -std=c++17 \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  "$SRC" -o "$BIN"

{
  date -u +%Y-%m-%dT%H:%M:%SZ
  git -C "$EXP" rev-parse HEAD
  shasum -a 256 "$SRC" "$BIN"
  printf '%s\n' \
    "$BIN $EXP/raw/driven_T2_T8 2 8 0.1 0.001 500 1000 0.01 64 2026090701 8" \
    "$BIN $EXP/raw/equilibrium_T5_T5 5 5 0.1 0.001 500 1000 0.01 64 2026090702 8" \
    "python3 $EXP/analyze_autocorr.py --experiment $EXP"
} > "$EXP/provenance/production_commands.txt"

mkdir -p "$EXP/raw/driven_T2_T8" "$EXP/raw/equilibrium_T5_T5"

"$BIN" "$EXP/raw/driven_T2_T8" \
  2 8 0.1 0.001 500 1000 0.01 64 2026090701 8 \
  > "$EXP/raw/driven_T2_T8/run.log" 2>&1

"$BIN" "$EXP/raw/equilibrium_T5_T5" \
  5 5 0.1 0.001 500 1000 0.01 64 2026090702 8 \
  > "$EXP/raw/equilibrium_T5_T5/run.log" 2>&1

python3 "$EXP/analyze_autocorr.py" --experiment "$EXP" \
  > "$EXP/analysis/analysis.log" 2>&1

find "$EXP/raw/driven_T2_T8" "$EXP/raw/equilibrium_T5_T5" \
  -type f -print0 | sort -z | xargs -0 shasum -a 256 \
  > "$EXP/provenance/raw_manifest.sha256"

date -u +%Y-%m-%dT%H:%M:%SZ > "$EXP/provenance/completed_utc.txt"
