#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_total_kde_2026-09-02"
BINARY="$ROOT/experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq"
PREFIX="$EXP/raw/driven"
FIFO="${PREFIX}_blocks.csv"
PARTIAL="${PREFIX}_blocks.csv.partial.zst"
ARCHIVE="${PREFIX}_blocks.csv.zst"
LOG="$EXP/raw/driven.log"
EXPECTED_RAW_SHA="4f728b3d0e007d704d90734b0888c00ec05b60f09385c6cfd079f3417d7a088f"

if [[ "$(pmset -g batt | head -n 1)" != *"AC Power"* ]]; then
  echo "refusing production restore while not on AC power" >&2
  exit 3
fi
if [[ ! -x "$BINARY" ]]; then
  echo "missing binary: $BINARY" >&2
  exit 2
fi
if [[ -e "$ARCHIVE" || -e "$PARTIAL" ]]; then
  echo "refusing to overwrite existing driven archive" >&2
  exit 2
fi

mkdir -p "$EXP/raw" "$EXP/provenance"
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
  [[ -p "$FIFO" ]] && unlink "$FIFO"
  exit "$status"
}
trap cleanup EXIT INT TERM

zstd -T0 -12 -q -c "$FIFO" > "$PARTIAL" &
compressor_pid=$!
command=(
  "$BINARY" sample_n3 10 2 3 8 500 20 7813 0.0005
  2026083133 8 "$PREFIX" 1
)
{
  printf 'started_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'command:'
  printf ' %q' "${command[@]}"
  printf '\n'
} | tee "$LOG"
"${command[@]}" 2>&1 | tee -a "$LOG"
wait "$compressor_pid"
compressor_pid=""
mv "$PARTIAL" "$ARCHIVE"
unlink "$FIFO"
trap - EXIT INT TERM

zstd -t -q "$ARCHIVE"
actual_raw_sha="$(zstd -dc "$ARCHIVE" | shasum -a 256 | awk '{print $1}')"
if [[ "$actual_raw_sha" != "$EXPECTED_RAW_SHA" ]]; then
  echo "restored raw SHA mismatch: $actual_raw_sha != $EXPECTED_RAW_SHA" | tee -a "$LOG" >&2
  exit 4
fi
{
  printf 'completed_utc=%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'compressed_sha256=%s\n' "$(shasum -a 256 "$ARCHIVE" | awk '{print $1}')"
  printf 'decompressed_sha256=%s\n' "$actual_raw_sha"
  printf 'binary_sha256=%s\n' "$(shasum -a 256 "$BINARY" | awk '{print $1}')"
} | tee -a "$LOG"

