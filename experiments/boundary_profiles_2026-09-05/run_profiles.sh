#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BINARY="$ROOT/bin/NLS_boundary_profiles"
MATRIX="$ROOT/RUN_MATRIX.csv"
STAGE="${1:-all}"

if [[ "$STAGE" != "gate" && "$STAGE" != "production" && "$STAGE" != "all" ]]; then
  echo "usage: $0 [gate|production|all]" >&2
  exit 2
fi

wait_for_ac() {
  while ! pmset -g batt | head -1 | grep -q "AC Power"; do
    echo "$(date -u +%FT%TZ) waiting_for_ac"
    sleep 60
  done
}

run_one() {
  local stage="$1" bc_label="$2" bc_id="$3" temp_label="$4"
  local T1="$5" Tn="$6" n="$7" replicate="$8" batches="$9"
  local burnin="${10}" measure="${11}" seed="${12}"
  local prefix="$ROOT/raw/${bc_label}_${temp_label}_n${n}_rep${replicate}"
  local log="$ROOT/logs/${bc_label}_${temp_label}_n${n}_rep${replicate}.log"
  local summary="${prefix}_summary.csv"

  if [[ -s "$summary" && -s "${prefix}_profile.csv" ]]; then
    echo "$(date -u +%FT%TZ) skip_complete prefix=$prefix"
    return 0
  fi
  wait_for_ac
  echo "$(date -u +%FT%TZ) start prefix=$prefix"
  echo "$BINARY $bc_id $T1 $Tn $n $batches $burnin $measure $seed $prefix"
  "$BINARY" "$bc_id" "$T1" "$Tn" "$n" "$batches" "$burnin" \
    "$measure" "$seed" "$prefix" > "$log" 2>&1
  test -s "$summary"
  test -s "${prefix}_profile.csv"
  echo "$(date -u +%FT%TZ) complete prefix=$prefix"
}

run_stage() {
  local wanted="$1"
  while IFS=, read -r stage bc_label bc_id temp_label T1 Tn n replicate \
      batches burnin measure seed; do
    [[ "$stage" == "stage" ]] && continue
    [[ "$stage" == "$wanted" ]] || continue
    run_one "$stage" "$bc_label" "$bc_id" "$temp_label" "$T1" "$Tn" \
      "$n" "$replicate" "$batches" "$burnin" "$measure" "$seed"
    if [[ "$wanted" == "gate" && "$replicate" == "1" ]]; then
      python3 "$ROOT/check_bc1_gate.py" --n "$n"
    fi
  done < "$MATRIX"
}

if [[ "$STAGE" == "gate" || "$STAGE" == "all" ]]; then
  run_stage gate
fi

if [[ "$STAGE" == "production" || "$STAGE" == "all" ]]; then
  for n in 25 50 100; do
    test "$(python3 -c 'import json,sys; print(str(json.load(open(sys.argv[1]))["pass"]).lower())' \
      "$ROOT/analysis/bc1_gate_n${n}.json")" = true
  done
  run_stage production
fi
