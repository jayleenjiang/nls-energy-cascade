#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_neep_data_2026-09-03"
mkdir -p "$EXP/analysis"

python3 "$EXP/analyze_neep_data.py" \
  --raw "$EXP/raw/driven_blocks.csv.zst" \
  --neep "$EXP/raw/driven_neep_transitions.csv.zst" \
  --label driven --tl 10 --tr 2 --output "$EXP/analysis" \
  > "$EXP/analysis/driven_analysis.log"

python3 "$EXP/analyze_neep_data.py" \
  --raw "$EXP/raw/equilibrium_blocks.csv.zst" \
  --neep "$EXP/raw/equilibrium_neep_transitions.csv.zst" \
  --label equilibrium --tl 6 --tr 6 --output "$EXP/analysis" \
  > "$EXP/analysis/equilibrium_analysis.log"
