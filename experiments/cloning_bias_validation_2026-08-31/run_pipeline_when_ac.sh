#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
root="$repo_root/experiments/cloning_bias_validation_2026-08-31"
lock="$root/pipeline.lock"

if ! mkdir "$lock" 2>/dev/null; then
  echo "Pipeline lock exists; refusing duplicate launch: $lock" >&2
  exit 1
fi
trap 'rmdir "$lock" 2>/dev/null || true' EXIT

if pgrep -f "$root/bin/entropy_cloning_v2 (controlled|endpoints)" >/dev/null 2>&1; then
  echo "Matching simulator process already active; refusing duplicate launch." >&2
  exit 1
fi

echo $$ > "$root/pipeline.pid"
date -u +%FT%TZ > "$root/STARTED_UTC.txt"

bash "$root/run_validation.sh" all-primary
python3 "$root/analyze_bias_validation.py" --stage plateau

python3 - "$root/analysis/optional_population_requests_primary.csv" <<'PY' > "$root/analysis/optional_targets.txt"
import csv, sys
with open(sys.argv[1], newline="") as handle:
    for row in csv.DictReader(handle):
        if int(row["run_8192"]):
            print(row["study"], row["n"])
PY

while read -r study n; do
  [[ -n "$study" ]] || continue
  if [[ "$study" == "n2_known_answer" ]]; then
    bash "$root/run_validation.sh" n2-8192
  else
    bash "$root/run_validation.sh" long-8192 "$n"
  fi
done < "$root/analysis/optional_targets.txt"

python3 "$root/analyze_bias_validation.py" --stage final
date -u +%FT%TZ > "$root/COMPLETED_UTC.txt"
