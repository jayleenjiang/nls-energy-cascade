#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_total_kde_2026-09-02"
DRIVEN="$EXP/raw/driven_blocks.csv.zst"
T6="$ROOT/experiments/entropy_ft_n3_equilibrium_2026-09-01/raw/T6_blocks.csv.zst"
T10="$ROOT/experiments/entropy_ft_n3_equilibrium_2026-09-01/raw/T10_blocks.csv.zst"
OUTPUT="$EXP/analysis"
LOG="$EXP/analysis.log"

for path in "$DRIVEN" "$T6" "$T10"; do
  if [[ ! -f "$path" ]]; then
    echo "missing input archive: $path" >&2
    exit 2
  fi
done
if [[ -e "$OUTPUT/FINAL_AUDIT.json" ]]; then
  echo "refusing to overwrite completed analysis" >&2
  exit 2
fi
mkdir -p "$OUTPUT"

command=(
  python3 "$EXP/analyze_kde_ft.py"
  --driven "$DRIVEN"
  --equilibrium-t6 "$T6"
  --equilibrium-t10 "$T10"
  --output "$OUTPUT"
  --write-derived
)
{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'command:'
  printf ' %q' "${command[@]}"
  printf '\n'
} | tee "$LOG"
"${command[@]}" 2>&1 | tee -a "$LOG"
printf 'completed_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "$LOG"

