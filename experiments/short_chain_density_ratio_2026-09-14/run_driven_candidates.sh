#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
PY=/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11
CACHE=/Users/jayleenjiang/NLS_section4_density_cache/short_chain_section4_2026-09-14
cd "$ROOT"
mkdir -p driven_fine_candidates provenance logs

architectures=$("$PY" - <<'PY'
import json
for value in json.load(open('EQUILIBRIUM_ADMISSIBLE_SET.json'))['admissible_architectures']:
    print(value)
PY
)

for architecture in $architectures; do
  for seed in 5201 5202 5203; do
    output="driven_fine_candidates/${architecture}_seed${seed}"
    if [[ -s "$output/validation_metrics.json" ]]; then
      continue
    fi
    command=("$PY" scripts/train_ratio.py --phase driven --timestep fine
      --architecture "$architecture" --seed "$seed" --output "$output")
    printf 'NLS_SECTION4_RAW_CACHE=%q ' "$CACHE" >> provenance/DRIVEN_FINE_COMMANDS.txt
    printf '%q ' "${command[@]}" >> provenance/DRIVEN_FINE_COMMANDS.txt
    printf '\n' >> provenance/DRIVEN_FINE_COMMANDS.txt
    NLS_SECTION4_RAW_CACHE="$CACHE" "${command[@]}" \
      > "logs/driven_fine_${architecture}_seed${seed}.log" 2>&1
  done
done
