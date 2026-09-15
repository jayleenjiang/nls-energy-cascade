#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BIN="$ROOT/bin/NLS_flux_canonical_profile"
LOG="$ROOT/pipeline.log"

verify_hashes() {
  local source_hash binary_hash
  source_hash="$(shasum -a 256 "$ROOT/src/NLS_flux_canonical_profile.cpp" | awk '{print $1}')"
  binary_hash="$(shasum -a 256 "$BIN" | awk '{print $1}')"
  [[ "$source_hash" == "c2b897c56a31b0ac38477df31a0900d58d68262f807eba5a38b6cf73e4eac0f5" ]]
  [[ "$binary_hash" == "e950c8b60a2b42e694bed0f19379d485f5b4e33ad160738e3f49fa0c6a944629" ]]
}

on_ac() {
  pmset -g batt | grep -q "AC Power"
}

wait_for_ac() {
  while ! on_ac; do
    printf '%s waiting for AC power\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG"
    sleep 300
  done
}

run_one() {
  local n="$1"
  local burnin="$2"
  local seed="$3"
  local rep="$4"
  local prefix="$ROOT/raw/canonical_T6_T6_n${n}_rep${rep}"
  if [[ -s "${prefix}_summary.csv" && -s "${prefix}_profile.csv" && -s "${prefix}_trajectory_profiles.csv" ]]; then
    printf '%s skip complete %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$prefix" >> "$LOG"
    return 0
  fi
  wait_for_ac
  printf '%s start n=%s rep=%s seed=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$n" "$rep" "$seed" >> "$LOG"
  "$BIN" 6 6 "$n" 16 "$burnin" 2000 0.0005 "$seed" 4 "$prefix" >> "$LOG" 2>&1
  printf '%s finish n=%s rep=%s seed=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$n" "$rep" "$seed" >> "$LOG"
}

verify_hashes
wait_for_ac
run_one 25 2000 2026090625 0 &
p1=$!
run_one 25 2000 2026090626 1 &
p2=$!
wait "$p1"
wait "$p2"

wait_for_ac
run_one 100 32000 2026090700 0 &
p3=$!
run_one 100 32000 2026090701 1 &
p4=$!
wait "$p3"
wait "$p4"

printf '%s all production complete\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG"
python3 "$ROOT/analyze.py" >> "$LOG" 2>&1
printf '0\n' > "$ROOT/pipeline.exitcode"
