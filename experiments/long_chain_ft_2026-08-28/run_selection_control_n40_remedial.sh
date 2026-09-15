#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
binary="$repo_root/flux/entropy_cloning_v2"
root="$repo_root/experiments/long_chain_ft_2026-08-28/numerical_controls_n40_remedial/N512_dt5e4_selection1"

if [[ ! -x "$binary" ]]; then
  echo "Build $binary first." >&2
  exit 1
fi

mkdir -p "$root"

run_one() {
  local k=$1 seed=$2 prefix=$3
  if [[ -s "${prefix}_summary.csv" && -s "${prefix}_timeseries.csv" ]]; then
    echo "skip completed: $prefix"
    return
  fi
  "$binary" controlled 10 2 40 512 500 120 1 0.0005 \
    "$k" 0.1 0.5 "$seed" 5 "$prefix"
}

for run in 1 2 3 4; do
  run_one 0.4 $((90300 + 2 * run - 1)) \
    "$root/n40_t120_k0p4_run${run}" &
  low_pid=$!
  run_one 0.6 $((90300 + 2 * run)) \
    "$root/n40_t120_k0p6_run${run}" &
  high_pid=$!
  wait "$low_pid" "$high_pid"
done
