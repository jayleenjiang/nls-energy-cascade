#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
binary="$repo_root/flux/entropy_cloning_v2"
root="$repo_root/experiments/long_chain_ft_phase2_2026-08-29/raw"
script_path="$repo_root/experiments/long_chain_ft_phase2_2026-08-29/run_phase2.sh"

if [[ ! -x "$binary" ]]; then
  echo "Build and self-test $binary first." >&2
  exit 1
fi

run_one() {
  local n=$1 clones=$2 horizon=$3 selection=$4 dt=$5
  local k=$6 seed=$7 prefix=$8
  if [[ -s "${prefix}_summary.csv" && -s "${prefix}_timeseries.csv" ]]; then
    echo "already complete: $prefix"
    return
  fi
  mkdir -p "$(dirname "$prefix")"
  {
    /usr/bin/time -p "$binary" controlled 10 2 "$n" "$clones" 500 \
      "$horizon" "$selection" "$dt" "$k" 0.1 0.5 "$seed" 5 "$prefix"
  } >"${prefix}.log" 2>&1
}

run_matrix() {
  local label=$1 n=$2 clones=$3 horizon=$4 selection=$5 dt=$6 seed_base=$7
  local output="$root/$label"
  for run in 1 2 3 4; do
    run_one "$n" "$clones" "$horizon" "$selection" "$dt" 0.4 \
      $((seed_base + 2 * run - 1)) "$output/n${n}_k0p4_run${run}" &
    low_pid=$!
    run_one "$n" "$clones" "$horizon" "$selection" "$dt" 0.6 \
      $((seed_base + 2 * run)) "$output/n${n}_k0p6_run${run}" &
    high_pid=$!
    set +e
    wait "$low_pid"
    low_status=$?
    wait "$high_pid"
    high_status=$?
    set -e
    if [[ $low_status -ne 0 || $high_status -ne 0 ]]; then
      echo "Phase-II pair failed in $label run $run; no automatic retry." >&2
      exit 1
    fi
  done
}

target=${1:-}
case "$target" in
  n20-dt)
    run_matrix n20_dt_N1024 20 1024 60 2 0.00025 91000
    ;;
  n20-selection)
    run_matrix n20_selection1_N1024 20 1024 60 1 0.0005 91100
    ;;
  n30-population)
    run_matrix n30_population_N4096 30 4096 60 2 0.0005 92000
    ;;
  n30-dt)
    run_matrix n30_dt_N4096 30 4096 60 2 0.00025 92100
    ;;
  n30-selection)
    run_matrix n30_selection1_N4096 30 4096 60 1 0.0005 92200
    ;;
  n40-selection)
    run_matrix n40_selection1_N1024 40 1024 120 1 0.0005 93000
    ;;
  n40-dt)
    run_matrix n40_dt_N1024 40 1024 120 2 0.00025 93100
    ;;
  stage1)
    "$script_path" n20-dt
    "$script_path" n20-selection
    "$script_path" n30-population
    "$script_path" n40-selection
    "$script_path" n40-dt
    ;;
  n30-controls)
    "$script_path" n30-dt
    "$script_path" n30-selection
    ;;
  *)
    echo "usage: $0 {stage1|n20-dt|n20-selection|n30-population|n30-controls|n30-dt|n30-selection|n40-selection|n40-dt}" >&2
    exit 2
    ;;
esac
