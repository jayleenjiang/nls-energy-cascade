#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "usage: $0 T seed label" >&2
  exit 2
fi

T="$1"
SEED="$2"
LABEL="$3"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_equilibrium_2026-09-01"
BINARY="$EXP/bin/entropy_ft_n3_eq"
PREFIX="$EXP/raw/${LABEL}"
FIFO="${PREFIX}_blocks.csv"
PARTIAL="${PREFIX}_blocks.csv.partial.zst"
ARCHIVE="${PREFIX}_blocks.csv.zst"
LOG="$EXP/raw/${LABEL}.log"

if [[ ! -x "$BINARY" ]]; then
  echo "missing binary: $BINARY" >&2
  exit 2
fi
if [[ -e "$ARCHIVE" || -e "$PARTIAL" ]]; then
  echo "refusing to overwrite existing output for $LABEL" >&2
  exit 2
fi

mkdir -p "$EXP/raw"
if [[ -p "$FIFO" ]]; then
  unlink "$FIFO"
elif [[ -e "$FIFO" ]]; then
  echo "refusing to replace non-FIFO path: $FIFO" >&2
  exit 2
fi
mkfifo "$FIFO"

compressor_pid=""
cleanup() {
  status=$?
  if [[ -n "$compressor_pid" ]] && kill -0 "$compressor_pid" 2>/dev/null; then
    kill "$compressor_pid" 2>/dev/null || true
    wait "$compressor_pid" 2>/dev/null || true
  fi
  if [[ -p "$FIFO" ]]; then
    unlink "$FIFO"
  fi
  exit "$status"
}
trap cleanup EXIT INT TERM

zstd -T0 -12 -q -c "$FIFO" > "$PARTIAL" &
compressor_pid=$!

command=(
  "$BINARY" sample_n3 "$T" "$T" 3 8 500 20 7813 0.0005
  "$SEED" 8 "$PREFIX" 1
)
printf 'command:' | tee "$LOG"
printf ' %q' "${command[@]}" | tee -a "$LOG"
printf '\n' | tee -a "$LOG"
"${command[@]}" 2>&1 | tee -a "$LOG"
wait "$compressor_pid"
compressor_pid=""
mv "$PARTIAL" "$ARCHIVE"
unlink "$FIFO"
trap - EXIT INT TERM

zstd -t -q "$ARCHIVE"
echo "complete: $ARCHIVE" | tee -a "$LOG"
