#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PYTHON=/Library/Frameworks/Python.framework/Versions/3.11/bin/python3.11
OUTROOT="$ROOT/stationary_density/equilibrium_candidates"
PROV="$ROOT/provenance"
mkdir -p "$OUTROOT" "$PROV"

export TF_CPP_MIN_LOG_LEVEL=1
export TF_NUM_INTRAOP_THREADS=4
export TF_NUM_INTEROP_THREADS=1
export OMP_NUM_THREADS=4
export OPENBLAS_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=4

date -u +%Y-%m-%dT%H:%M:%SZ > "$PROV/equilibrium_candidates_started_utc.txt"
: > "$PROV/equilibrium_candidate_commands.txt"

for architecture in "8 64,64" "16 64,64" "16 128,128"; do
  read -r components widths <<< "$architecture"
  width_tag="${widths/,/x}"
  for lambda_fp in 0.01 0.1; do
    lambda_tag="${lambda_fp/./p}"
    for seed in 4101 4102 4103; do
      name="K$(printf '%02d' "$components")_w${width_tag}_l${lambda_tag}_seed${seed}"
      out="$OUTROOT/$name"
      if [[ -s "$out/validation_metrics.json" ]]; then
        echo "skip completed $name"
        continue
      fi
      mkdir -p "$out"
      command=(
        "$PYTHON" "$ROOT/scripts/train_density.py"
        --case equilibrium_dt2p5e-4
        --components "$components"
        --widths "$widths"
        --lambda-fp "$lambda_fp"
        --seed "$seed"
        --output "$out"
        --pretrain-epochs 50
        --joint-epochs 450
        --steps-per-epoch 128
        --nll-batch 4096
        --fp-batch 512
        --gmm-sample 200000
        --validation-sample 200000
        --validation-fp-sample 20000
      )
      printf '%q ' "${command[@]}" | tee -a "$PROV/equilibrium_candidate_commands.txt"
      printf '\n' | tee -a "$PROV/equilibrium_candidate_commands.txt"
      "${command[@]}" 2>&1 | tee "$out/run.log"
    done
  done
done

"$PYTHON" "$ROOT/scripts/select_density_candidate.py" \
  2>&1 | tee "$OUTROOT/selection.log"
date -u +%Y-%m-%dT%H:%M:%SZ > "$PROV/equilibrium_candidates_completed_utc.txt"

