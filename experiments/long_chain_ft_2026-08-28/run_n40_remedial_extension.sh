#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
binary="$repo_root/flux/entropy_cloning_v2"
root="$repo_root/experiments/long_chain_ft_2026-08-28/selection2_n40_k04_06_remedial"
script_path="$repo_root/experiments/long_chain_ft_2026-08-28/run_n40_remedial_extension.sh"

if [[ ! -x "$binary" ]]; then
  echo "Build $binary first." >&2
  exit 1
fi

run_one() {
  local population=$1 k=$2 seed=$3 prefix=$4
  if [[ -s "${prefix}_summary.csv" && -s "${prefix}_timeseries.csv" ]]; then
    echo "skip completed: $prefix"
    return
  fi
  "$binary" controlled 10 2 40 "$population" 500 120 2 0.0005 \
    "$k" 0.1 0.5 "$seed" 5 "$prefix"
}

run_matrix() {
  local population=$1
  local seed_base=$2
  local output="$root/N${population}_t120"
  mkdir -p "$output"
  for run in 1 2 3 4; do
    run_one "$population" 0.4 $((seed_base + 2 * run - 1)) \
      "$output/n40_t120_k0p4_run${run}" &
    local low_pid=$!
    run_one "$population" 0.6 $((seed_base + 2 * run)) \
      "$output/n40_t120_k0p6_run${run}" &
    local high_pid=$!
    wait "$low_pid" "$high_pid"
  done
}

target=${1:-all}
case "$target" in
  N512)
    run_matrix 512 89600
    ;;
  N1024)
    run_matrix 1024 89700
    ;;
  all)
    "$script_path" N512
    "$script_path" N1024
    ;;
  *)
    echo "usage: $0 {N512|N1024|all}" >&2
    exit 2
    ;;
esac
