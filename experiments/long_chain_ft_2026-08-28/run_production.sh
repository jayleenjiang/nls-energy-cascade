#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
binary="$repo_root/flux/entropy_cloning_v2"
experiment_root="$repo_root/experiments/long_chain_ft_2026-08-28"
script_path="$experiment_root/run_production.sh"

if [[ ! -x "$binary" ]]; then
  echo "Build $binary from flux/NLS_entropy_cloning.cpp first." >&2
  exit 1
fi

run_one() {
  local n=$1 population=$2 horizon=$3 k=$4 seed=$5 threads=$6 prefix=$7
  if [[ -s "${prefix}_summary.csv" && -s "${prefix}_timeseries.csv" ]]; then
    echo "skip completed: $prefix"
    return
  fi
  "$binary" controlled 10 2 "$n" "$population" 500 "$horizon" \
    2 0.0005 "$k" 0.1 0.5 "$seed" "$threads" "$prefix"
}

run_pair_matrix() {
  local n=$1 population=$2 horizon=$3 k_low=$4 k_high=$5
  local seed_base=$6 threads=$7 output_dir=$8
  mkdir -p "$output_dir"
  for run in 1 2 3 4; do
    local seed_low=$((seed_base + 2 * run - 1))
    local seed_high=$((seed_base + 2 * run))
    local low_tag=${k_low/./p}
    local high_tag=${k_high/./p}
    run_one "$n" "$population" "$horizon" "$k_low" "$seed_low" \
      "$threads" "$output_dir/n${n}_t${horizon}_k${low_tag}_run${run}" &
    local low_pid=$!
    run_one "$n" "$population" "$horizon" "$k_high" "$seed_high" \
      "$threads" "$output_dir/n${n}_t${horizon}_k${high_tag}_run${run}" &
    local high_pid=$!
    wait "$low_pid" "$high_pid"
  done
}

target=${1:-all}
case "$target" in
  n10)
    run_pair_matrix 10 1024 80 0.3 0.7 87700 5 \
      "$experiment_root/selection2_n10/N1024_t80"
    run_pair_matrix 10 2048 60 0.3 0.7 87800 5 \
      "$experiment_root/selection2_n10/N2048_t60"
    ;;
  n20)
    run_pair_matrix 20 512 60 0.4 0.6 88300 5 \
      "$experiment_root/selection2_n20_k04_06/N512_t60"
    run_pair_matrix 20 1024 60 0.4 0.6 88400 5 \
      "$experiment_root/selection2_n20_k04_06/N1024_t60"
    ;;
  n30)
    run_pair_matrix 30 512 60 0.4 0.6 88500 3 \
      "$experiment_root/selection2_n30_k04_06/N512_t60"
    run_pair_matrix 30 1024 60 0.4 0.6 88700 5 \
      "$experiment_root/selection2_n30_k04_06/N1024_t60"
    ;;
  n40)
    run_pair_matrix 40 512 60 0.4 0.6 88600 2 \
      "$experiment_root/selection2_n40_k04_06/N512_t60"
    run_pair_matrix 40 1024 80 0.4 0.6 88800 5 \
      "$experiment_root/selection2_n40_k04_06/N1024_t80"
    ;;
  all)
    "$script_path" n10
    "$script_path" n20
    "$script_path" n30
    "$script_path" n40
    ;;
  *)
    echo "usage: $0 {n10|n20|n30|n40|all}" >&2
    exit 2
    ;;
esac
