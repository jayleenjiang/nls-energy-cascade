#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "usage: $0 label" >&2
  exit 2
fi

LABEL="$1"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_neep_data_2026-09-03"
RAW="$EXP/raw/${LABEL}_blocks.csv.zst"
PARTIAL="$EXP/raw/${LABEL}_neep_transitions.csv.partial.zst"
OUTPUT="$EXP/raw/${LABEL}_neep_transitions.csv.zst"
LOG="$EXP/raw/${LABEL}_augment.log"

if [[ ! -f "$RAW" ]]; then
  echo "missing raw archive: $RAW" >&2
  exit 2
fi
if [[ -e "$OUTPUT" || -e "$PARTIAL" ]]; then
  echo "refusing to overwrite NEEP archive for $LABEL" >&2
  exit 2
fi

set -o pipefail
zstd -dc "$RAW" \
  | python3 "$EXP/augment_neep.py" 2> "$LOG" \
  | zstd -T0 -12 -q -c > "$PARTIAL"
mv "$PARTIAL" "$OUTPUT"
zstd -t -q "$OUTPUT"
echo "complete: $OUTPUT" | tee -a "$LOG"
