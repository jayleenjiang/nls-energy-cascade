#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
binary="$repo_root/flux/entropy_cloning_v2"
root="$repo_root/experiments/long_chain_ft_2026-08-28/numerical_controls_n10_remedial"
script_path="$repo_root/experiments/long_chain_ft_2026-08-28/run_numerical_controls_n10_remedial.sh"

if [[ ! -x "$binary" ]]; then
  echo "Build $binary first." >&2
  exit 1
fi

run_one() {
  local selection=$1 dt=$2 k=$3 seed=$4 prefix=$5
  if [[ -s "${prefix}_summary.csv" && -s "${prefix}_timeseries.csv" ]]; then
    echo "skip completed: $prefix"
    return
  fi
  "$binary" controlled 10 2 10 1024 500 60 "$selection" "$dt" \
    "$k" 0.1 0.5 "$seed" 5 "$prefix"
}

run_matrix() {
  local label=$1 selection=$2 dt=$3 seed_base=$4
  local output="$root/$label"
  mkdir -p "$output"
  for run in 1 2 3 4; do
    run_one "$selection" "$dt" 0.3 $((seed_base + 2 * run - 1)) \
      "$output/n10_t60_k0p3_run${run}" &
    low_pid=$!
    run_one "$selection" "$dt" 0.7 $((seed_base + 2 * run)) \
      "$output/n10_t60_k0p7_run${run}" &
    high_pid=$!
    wait "$low_pid" "$high_pid"
  done
}

target=${1:-all}
case "$target" in
  timestep)
    run_matrix N1024_dt2p5e4_selection2 2 0.00025 90100
    ;;
  selection)
    run_matrix N1024_dt5e4_selection1 1 0.0005 90200
    ;;
  all)
    "$script_path" timestep
    "$script_path" selection
    ;;
  *)
    echo "usage: $0 {timestep|selection|all}" >&2
    exit 2
    ;;
esac
