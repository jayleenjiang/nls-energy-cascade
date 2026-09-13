#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
phase3="$repo_root/experiments/long_chain_ft_multik_2026-08-31"
binary="$repo_root/flux/entropy_cloning_v2"
raw="$phase3/raw"

if [[ ! -x "$binary" ]]; then
  echo "Build and self-test $binary first." >&2
  exit 1
fi

k_tag() {
  printf '%s' "$1" | sed 's/^0\./0p/'
}

run_one() {
  local n=$1 clones=$2 horizon=$3 selection=$4 dt=$5
  local k=$6 seed=$7 prefix=$8
  if [[ -s "${prefix}_summary.csv" && -s "${prefix}_timeseries.csv" ]]; then
    echo "already complete: $prefix"
    return
  fi
  if [[ -e "${prefix}_summary.csv" || -e "${prefix}_timeseries.csv" ]]; then
    echo "partial output exists; refusing to overwrite: $prefix" >&2
    return 1
  fi
  mkdir -p "$(dirname "$prefix")"
  {
    /usr/bin/time -p "$binary" controlled 10 2 "$n" "$clones" 500 \
      "$horizon" "$selection" "$dt" "$k" 0.1 0.5 "$seed" 5 "$prefix"
  } >"${prefix}.log" 2>&1
}

run_pair() {
  local n=$1 clones=$2 horizon=$3 selection=$4 dt=$5
  local k_low=$6 k_high=$7 seed_low=$8 seed_high=$9 output=${10}
  local low_tag high_tag low_prefix high_prefix
  low_tag=$(k_tag "$k_low")
  high_tag=$(k_tag "$k_high")
  low_prefix="$output/n${n}_k${low_tag}_s${seed_low}"
  high_prefix="$output/n${n}_k${high_tag}_s${seed_high}"
  run_one "$n" "$clones" "$horizon" "$selection" "$dt" \
    "$k_low" "$seed_low" "$low_prefix" &
  local low_pid=$!
  run_one "$n" "$clones" "$horizon" "$selection" "$dt" \
    "$k_high" "$seed_high" "$high_prefix" &
  local high_pid=$!
  set +e
  wait "$low_pid"
  local low_status=$?
  wait "$high_pid"
  local high_status=$?
  set -e
  if [[ $low_status -ne 0 || $high_status -ne 0 ]]; then
    echo "pair failed: n=$n, k=($k_low,$k_high)" >&2
    return 1
  fi
}

run_pilot() {
  while IFS=, read -r n k_low k_high seed_low seed_high clones horizon selection dt; do
    [[ $n == n ]] && continue
    run_pair "$n" "$clones" "$horizon" "$selection" "$dt" \
      "$k_low" "$k_high" "$seed_low" "$seed_high" \
      "$raw/pilot/n${n}_k$(k_tag "$k_low")_$(k_tag "$k_high")"
  done < "$phase3/PILOT_MATRIX.csv"
}

run_production() {
  [[ -s "$phase3/FROZEN_GRID.csv" ]] || {
    echo "FROZEN_GRID.csv is required." >&2
    exit 1
  }
  while IFS=, read -r n role k_low k_high clones horizon baseline_low \
      baseline_high production_seed_base population_seed_base; do
    [[ $n == n ]] && continue
    local output="$raw/production/n${n}_${role}_k$(k_tag "$k_low")_$(k_tag "$k_high")"
    for run in 1 2 3 4; do
      run_pair "$n" "$clones" "$horizon" 2 0.0005 \
        "$k_low" "$k_high" \
        $((production_seed_base + 2 * run - 1)) \
        $((production_seed_base + 2 * run)) "$output"
    done
  done < "$phase3/FROZEN_GRID.csv"
}

run_population_controls() {
  [[ -s "$phase3/FROZEN_GRID.csv" ]] || exit 1
  while IFS=, read -r n role k_low k_high clones horizon baseline_low \
      baseline_high production_seed_base population_seed_base; do
    [[ $n == n || $role != outer ]] && continue
    local lower_clones=$((clones / 2))
    local output="$raw/population/n${n}_${role}_k$(k_tag "$k_low")_$(k_tag "$k_high")_N${lower_clones}"
    for run in 1 2 3 4; do
      run_pair "$n" "$lower_clones" "$horizon" 2 0.0005 \
        "$k_low" "$k_high" \
        $((population_seed_base + 2 * run - 1)) \
        $((population_seed_base + 2 * run)) "$output"
    done
  done < "$phase3/FROZEN_GRID.csv"
}

run_n40_controls() {
  local row
  row=$(awk -F, '$1==40 && $2=="outer" {print; exit}' "$phase3/FROZEN_GRID.csv")
  [[ -n $row ]] || {
    echo "missing n=40 outer pair" >&2
    exit 1
  }
  IFS=, read -r n role k_low k_high clones horizon baseline_low \
    baseline_high production_seed_base population_seed_base <<< "$row"
  local dt_output="$raw/n40_controls/timestep_k$(k_tag "$k_low")_$(k_tag "$k_high")"
  local sel_output="$raw/n40_controls/selection_k$(k_tag "$k_low")_$(k_tag "$k_high")"
  for run in 1 2 3 4; do
    run_pair 40 "$clones" "$horizon" 2 0.00025 \
      "$k_low" "$k_high" $((984000 + 2 * run - 1)) \
      $((984000 + 2 * run)) "$dt_output"
    run_pair 40 "$clones" "$horizon" 1 0.0005 \
      "$k_low" "$k_high" $((994000 + 2 * run - 1)) \
      $((994000 + 2 * run)) "$sel_output"
  done
}

target=${1:-}
case "$target" in
  pilot) run_pilot ;;
  production) run_production ;;
  population) run_population_controls ;;
  n40-controls) run_n40_controls ;;
  *)
    echo "usage: $0 {pilot|production|population|n40-controls}" >&2
    exit 2
    ;;
esac

