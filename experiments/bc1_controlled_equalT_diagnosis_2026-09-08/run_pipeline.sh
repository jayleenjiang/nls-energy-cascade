#!/bin/bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")" && pwd)"
BIN="$ROOT/bin/NLS_bc1_cartesian_profile"
SRC="$ROOT/src/NLS_bc1_cartesian_profile.cpp"
HASHFILE="$ROOT/FROZEN_HASHES.sha256"
COMMANDS="$ROOT/provenance/COMMANDS.csv"
MAX_PARALLEL=8

mkdir -p "$ROOT/raw" "$ROOT/analysis" "$ROOT/provenance"
cd "$ROOT"

if [[ ! -x "$BIN" ]]; then
  echo "missing production binary: $BIN" >&2
  exit 2
fi
shasum -a 256 -c "$HASHFILE"
git -C "$ROOT" rev-parse HEAD > "$ROOT/provenance/LAUNCH_COMMIT.txt"
date -u +%Y-%m-%dT%H:%M:%SZ > "$ROOT/provenance/STARTED_UTC.txt"
if [[ ! -e "$COMMANDS" ]]; then
  echo 'n,base_seed,stream_offset,command' > "$COMMANDS"
fi

run_length() {
  local n="$1" burn="$2" seed_a="$3" seed_b="$4"
  local -a specs=()
  local seed offset
  for seed in "$seed_a" "$seed_b"; do
    for offset in $(seq 0 16 240); do
      specs+=("$seed:$offset")
    done
  done

  local index=0
  while (( index < ${#specs[@]} )); do
    local -a pids=()
    local -a labels=()
    local slot
    for ((slot=0; slot<MAX_PARALLEL && index<${#specs[@]}; slot++,index++)); do
      IFS=: read -r seed offset <<< "${specs[$index]}"
      local prefix="$ROOT/raw/n${n}_seed${seed}_offset$(printf '%03d' "$offset")"
      if [[ -s "${prefix}_summary.csv" ]]; then
        echo "skip completed n=$n seed=$seed offset=$offset"
        continue
      fi
      if compgen -G "${prefix}_*" > /dev/null; then
        echo "refusing to overwrite incomplete output for $prefix" >&2
        exit 3
      fi
      local command="$BIN $prefix 6 6 $n 0.1 0.0005 $burn 2000 0.01 $seed $offset"
      printf '%s,%s,%s,"%s"\n' "$n" "$seed" "$offset" "$command" >> "$COMMANDS"
      echo "launch n=$n seed=$seed offset=$offset"
      "$BIN" "$prefix" 6 6 "$n" 0.1 0.0005 "$burn" 2000 0.01 \
        "$seed" "$offset" > "${prefix}_run.log" 2>&1 &
      pids+=("$!")
      labels+=("n=$n seed=$seed offset=$offset")
    done

    local failed=0
    local i
    for i in "${!pids[@]}"; do
      if ! wait "${pids[$i]}"; then
        echo "FAILED ${labels[$i]}" >&2
        failed=1
      else
        echo "complete ${labels[$i]}"
      fi
    done
    if (( failed )); then
      python3 "$ROOT/analyze.py" --n "$n" || true
      exit 4
    fi
  done
}

run_length 25 2000 2026090825 2026090826
python3 "$ROOT/analyze.py" --n 25 --require-complete \
  > "$ROOT/analysis/n25_analysis_stdout.json"

python3 - "$ROOT/analysis/n25_summary.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]))
if not (x['porting']['pass'] and x['numerical_pass']):
    raise SystemExit('n=25 porting/numerical gate failed; n=100 is forbidden')
PY

run_length 100 32000 2026090900 2026090901
python3 "$ROOT/analyze.py" --n 100 --require-complete \
  > "$ROOT/analysis/n100_analysis_stdout.json"
python3 "$ROOT/analyze.py" --n all --require-complete \
  > "$ROOT/analysis/final_analysis_stdout.json"

find "$ROOT/raw" "$ROOT/analysis" -type f -print0 \
  | sort -z | xargs -0 shasum -a 256 > "$ROOT/provenance/OUTPUT_HASHES.sha256"
date -u +%Y-%m-%dT%H:%M:%SZ > "$ROOT/provenance/COMPLETED_UTC.txt"

