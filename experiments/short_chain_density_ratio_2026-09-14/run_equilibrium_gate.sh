#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
PY=/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11
cd "$ROOT"
mkdir -p equilibrium_candidates provenance logs

for architecture in linear mlp32 mlp64; do
  for seed in 5201 5202 5203; do
    output="equilibrium_candidates/${architecture}_seed${seed}"
    if [[ -s "$output/validation_metrics.json" ]]; then
      continue
    fi
    command=("$PY" scripts/train_ratio.py --phase equilibrium --timestep fine
      --architecture "$architecture" --seed "$seed" --output "$output")
    printf '%q ' "${command[@]}" >> provenance/EQUILIBRIUM_COMMANDS.txt
    printf '\n' >> provenance/EQUILIBRIUM_COMMANDS.txt
    "${command[@]}" > "logs/${architecture}_seed${seed}.log" 2>&1
  done
done

"$PY" scripts/select_equilibrium.py > logs/equilibrium_selection.log 2>&1

