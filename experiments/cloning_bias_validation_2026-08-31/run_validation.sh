#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
root="$repo_root/experiments/cloning_bias_validation_2026-08-31"
binary="$root/bin/entropy_cloning_v2"
matrix="$root/RUN_MATRIX.csv"

if [[ ! -x "$binary" ]]; then
  echo "Missing release binary: run build_and_selftest.sh first." >&2
  exit 1
fi
if [[ ! -s "$matrix" ]]; then
  echo "Missing RUN_MATRIX.csv: run generate_run_matrix.py first." >&2
  exit 1
fi

power_is_ac() {
  pmset -g batt | head -n 1 | grep -q "AC Power"
}

wait_for_ac() {
  while ! power_is_ac; do
    echo "$(date -u +%FT%TZ) waiting for AC power"
    sleep 300
  done
}

lookup_row() {
  local study=$1 stage=$2 n=$3 clones=$4 k=$5 run=$6
  python3 - "$matrix" "$study" "$stage" "$n" "$clones" "$k" "$run" <<'PY'
import csv, sys
path, study, stage, n, clones, k, run = sys.argv[1:]
matches = []
with open(path, newline="") as handle:
    for row in csv.DictReader(handle):
        if (row["study"] == study and row["stage"] == stage
                and int(row["n"]) == int(n)
                and int(row["clone_count"]) == int(clones)
                and abs(float(row["k"]) - float(k)) < 1e-12
                and int(row["run"]) == int(run)):
            matches.append(row)
if len(matches) != 1:
    raise SystemExit(f"expected one matrix row, found {len(matches)}")
r = matches[0]
print("\t".join(r[key] for key in (
    "n", "clone_count", "burnin", "horizon", "selection_time", "dt",
    "k", "gauge_shift", "control_scale", "seed", "threads", "prefix")))
PY
}

run_one() {
  local study=$1 stage=$2 n=$3 clones=$4 k=$5 run=$6
  local row
  row=$(lookup_row "$study" "$stage" "$n" "$clones" "$k" "$run")
  local burnin horizon selection dt gauge control seed threads relative
  IFS=$'\t' read -r n clones burnin horizon selection dt k gauge control seed threads relative <<<"$row"
  local prefix="$root/$relative"
  local summary="${prefix}_summary.csv"
  local timeseries="${prefix}_timeseries.csv"
  local log="${prefix}.log"
  local command_file="${prefix}.command.txt"
  mkdir -p "$(dirname "$prefix")"

  if [[ -s "$summary" && -s "$timeseries" ]]; then
    echo "already complete: $prefix"
    return 0
  fi
  if [[ -e "$summary" || -e "$timeseries" || -e "$log" ]]; then
    echo "Incomplete prior artifact exists; no automatic retry: $prefix" >&2
    return 1
  fi

  local cmd=("$binary" controlled 10 2 "$n" "$clones" "$burnin"
    "$horizon" "$selection" "$dt" "$k" "$gauge" "$control" "$seed"
    "$threads" "$prefix")
  printf '%q ' "${cmd[@]}" > "$command_file"
  printf '\n' >> "$command_file"
  {
    /usr/bin/time -p "${cmd[@]}"
  } >"$log" 2>&1
}

wait_pair() {
  local pid_a=$1 pid_b=$2 label=$3
  set +e
  wait "$pid_a"; local status_a=$?
  wait "$pid_b"; local status_b=$?
  set -e
  if [[ $status_a -ne 0 || $status_b -ne 0 ]]; then
    echo "Run pair failed ($label); no automatic retry." >&2
    exit 1
  fi
}

run_n2_reference() {
  wait_for_ac
  local prefix="$root/raw/n2_direct_reference/n2_t0p1"
  local summary="${prefix}_summary.csv"
  local blocks="${prefix}_blocks.csv"
  local log="${prefix}.log"
  local command_file="${prefix}.command.txt"
  mkdir -p "$(dirname "$prefix")"
  if [[ -s "$summary" && -s "$blocks" ]]; then
    echo "already complete: $prefix"
    return
  fi
  if [[ -e "$summary" || -e "$blocks" || -e "$log" ]]; then
    echo "Incomplete prior n=2 direct artifact exists; no automatic retry." >&2
    exit 1
  fi
  local cmd=("$binary" endpoints 10 2 2 1024 500 0.1 1024 0.0005
    210000001 10 "$prefix")
  printf '%q ' "${cmd[@]}" > "$command_file"
  printf '\n' >> "$command_file"
  {
    /usr/bin/time -p "${cmd[@]}"
  } >"$log" 2>&1
}

run_n2_population() {
  local stage=$1 clones=$2
  for run in 1 2 3 4; do
    wait_for_ac
    run_one n2_known_answer "$stage" 2 "$clones" 0.3 "$run" & local low_pid=$!
    run_one n2_known_answer "$stage" 2 "$clones" 0.7 "$run" & local high_pid=$!
    wait_pair "$low_pid" "$high_pid" "n2 N=$clones run=$run k=0.3/0.7"
  done
  for first in 1 3; do
    wait_for_ac
    run_one n2_known_answer "$stage" 2 "$clones" 0.5 "$first" & local pid_a=$!
    run_one n2_known_answer "$stage" 2 "$clones" 0.5 "$((first + 1))" & local pid_b=$!
    wait_pair "$pid_a" "$pid_b" "n2 N=$clones k=0.5 runs=$first,$((first + 1))"
  done
}

tilts_for_n() {
  case "$1" in
    10) echo "0.3 0.7" ;;
    20|30|40) echo "0.4 0.6" ;;
    *) echo "unsupported n=$1" >&2; exit 2 ;;
  esac
}

run_long_population() {
  local stage=$1 n=$2 clones=$3
  local low high
  read -r low high <<<"$(tilts_for_n "$n")"
  for run in 1 2 3 4; do
    wait_for_ac
    run_one long_chain "$stage" "$n" "$clones" "$low" "$run" & local low_pid=$!
    run_one long_chain "$stage" "$n" "$clones" "$high" "$run" & local high_pid=$!
    wait_pair "$low_pid" "$high_pid" "n=$n N=$clones run=$run"
  done
}

run_primary() {
  run_n2_reference
  for clones in 512 1024 2048 4096; do
    run_n2_population primary "$clones"
  done
  for n in 10 20 30 40; do
    for clones in 512 1024 2048 4096; do
      run_long_population primary "$n" "$clones"
    done
  done
}

target=${1:-}
case "$target" in
  n2-reference)
    run_n2_reference
    ;;
  n2-primary)
    for clones in 512 1024 2048 4096; do run_n2_population primary "$clones"; done
    ;;
  n2-8192)
    run_n2_population conditional_8192 8192
    ;;
  long-primary)
    for n in 10 20 30 40; do
      for clones in 512 1024 2048 4096; do
        run_long_population primary "$n" "$clones"
      done
    done
    ;;
  long-8192)
    [[ $# -eq 2 ]] || { echo "usage: $0 long-8192 n" >&2; exit 2; }
    run_long_population conditional_8192 "$2" 8192
    ;;
  all-primary)
    run_primary
    ;;
  *)
    echo "usage: $0 {n2-reference|n2-primary|n2-8192|long-primary|long-8192 n|all-primary}" >&2
    exit 2
    ;;
esac
