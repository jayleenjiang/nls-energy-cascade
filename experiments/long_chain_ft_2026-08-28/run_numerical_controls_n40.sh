#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
binary="$repo_root/flux/entropy_cloning_v2"
root="$repo_root/experiments/long_chain_ft_2026-08-28/numerical_controls_n40"
script_path="$repo_root/experiments/long_chain_ft_2026-08-28/run_numerical_controls_n40.sh"

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
  "$binary" controlled 10 2 40 512 500 120 "$selection" "$dt" \
    "$k" 0.1 0.5 "$seed" 5 "$prefix"
}

run_matrix() {
  local label=$1 selection=$2 dt=$3 seed_base=$4
  local output="$root/$label"
  mkdir -p "$output"
  for run in 1 2 3 4; do
    run_one "$selection" "$dt" 0.4 $((seed_base + 2 * run - 1)) \
      "$output/n40_t120_k0p4_run${run}" &
    local low_pid=$!
    run_one "$selection" "$dt" 0.6 $((seed_base + 2 * run)) \
      "$output/n40_t120_k0p6_run${run}" &
    local high_pid=$!
    wait "$low_pid" "$high_pid"
  done
}

target=${1:-all}
case "$target" in
  timestep)
    run_matrix N512_dt2p5e4_selection2 2 0.00025 89400
    ;;
  selection)
    run_matrix N512_dt5e4_selection4 4 0.0005 89500
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
