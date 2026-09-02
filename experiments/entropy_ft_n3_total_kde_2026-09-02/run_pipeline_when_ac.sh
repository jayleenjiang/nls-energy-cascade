#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_total_kde_2026-09-02"
LOG="$EXP/pipeline.log"
DRIVEN="$EXP/raw/driven_blocks.csv.zst"

mkdir -p "$EXP/raw"
echo "waiting_for_ac_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG"
while [[ "$(pmset -g batt | head -n 1)" != *"AC Power"* ]]; do
  sleep 60
done
echo "ac_detected_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG"

if [[ ! -f "$DRIVEN" ]]; then
  "$EXP/run_restore_driven.sh" >> "$LOG" 2>&1
else
  zstd -t -q "$DRIVEN"
  echo "using_existing_driven_archive=$DRIVEN" >> "$LOG"
fi

if [[ ! -f "$EXP/analysis/FINAL_AUDIT.json" ]]; then
  "$EXP/run_analysis.sh" >> "$LOG" 2>&1
else
  echo "analysis_already_complete" >> "$LOG"
fi
echo "pipeline_completed_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG"

