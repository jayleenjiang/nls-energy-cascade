#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
binary="$repo_root/flux/entropy_cloning_v2"
root="$repo_root/experiments/long_chain_ft_2026-08-28/selection2_n30_k04_06_remedial/N2048_t60"

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
  "$binary" controlled 10 2 30 2048 500 60 2 0.0005 \
    "$k" 0.1 0.5 "$seed" 5 "$prefix"
}

for run in 1 2 3 4; do
  run_one 0.4 $((89800 + 2 * run - 1)) \
    "$root/n30_t60_k0p4_run${run}" &
  low_pid=$!
  run_one 0.6 $((89800 + 2 * run)) \
    "$root/n30_t60_k0p6_run${run}" &
  high_pid=$!
  wait "$low_pid" "$high_pid"
done
