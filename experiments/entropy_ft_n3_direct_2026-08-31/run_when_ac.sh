#!/bin/bash
set -euo pipefail

EXP=/Users/jayleenjiang/Documents/NLS/experiments/entropy_ft_n3_direct_2026-08-31
LOG="$EXP/launcher.log"

mkdir -p "$EXP"
while ! pmset -g batt | grep -q "AC Power"; do
  printf '%s waiting-for-AC\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG"
  sleep 300
done

printf '%s AC-detected; starting frozen production\n' \
  "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG"
exec "$EXP/run_production.sh" >> "$LOG" 2>&1
