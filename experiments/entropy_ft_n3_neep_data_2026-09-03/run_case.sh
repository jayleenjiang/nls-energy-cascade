#!/usr/bin/env bash
set -euo pipefail

if [[ $# -ne 4 ]]; then
  echo "usage: $0 TL TR seed label" >&2
  exit 2
fi

TL="$1"
TR="$2"
SEED="$3"
LABEL="$4"
ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_neep_data_2026-09-03"
BINARY="$ROOT/experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq"
PREFIX="$EXP/raw/${LABEL}"
FIFO="${PREFIX}_blocks.csv"
PARTIAL="${PREFIX}_blocks.csv.partial.zst"
ARCHIVE="${PREFIX}_blocks.csv.zst"
LOG="$EXP/raw/${LABEL}.log"
EXPECTED_SOURCE_SHA="9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2"
EXPECTED_BINARY_SHA="93c4aa6d046f3c35cd0ea1136091fc44ef1b64baf6237e4663ca9e86b995a156"

if [[ ! -x "$BINARY" ]]; then
  echo "missing validated binary: $BINARY" >&2
  exit 2
fi
if [[ "$(shasum -a 256 "$ROOT/flux/NLS_entropy_ft.cpp" | awk '{print $1}')" != "$EXPECTED_SOURCE_SHA" ]]; then
  echo "validated source SHA mismatch" >&2
  exit 2
fi
if [[ "$(shasum -a 256 "$BINARY" | awk '{print $1}')" != "$EXPECTED_BINARY_SHA" ]]; then
  echo "validated binary SHA mismatch" >&2
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
  "$BINARY" sample_n3 "$TL" "$TR" 3 8 500 0.1 39063 0.0005
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
