#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SOURCE="$HERE/source/NLS_entropy_ft_1905cf.cpp"
BIN="$HERE/bin/entropy_ft_affinity_equilibrium"
RAW="$HERE/raw"
OUTPUT_PREFIX="$RAW/dbeta_0p000000_n10"
BLOCKS="$OUTPUT_PREFIX"_blocks.csv
LOG="$RAW/dbeta_0p000000_n10.log"
MANIFEST="$HERE/EQUILIBRIUM_PRODUCTION_MANIFEST.txt"
COMMANDS="$HERE/EQUILIBRIUM_COMMAND.tsv"
RAW_MANIFEST="$HERE/EQUILIBRIUM_RAW_DATA_MANIFEST.tsv"
CXX="${CXX:-clang++}"
EIGEN3_INCLUDE="${EIGEN3_INCLUDE:-/opt/homebrew/include/eigen3}"
LIBOMP_PREFIX="${LIBOMP_PREFIX:-/opt/homebrew/opt/libomp}"
SOURCE_COMMIT="1905cf4e606a4a7f4dd8930caa64bd4cc861e9d4"
SOURCE_SHA="98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3"

mkdir -p "$HERE/bin" "$RAW"

if ! pmset -g batt | head -1 | grep -q "AC Power"; then
  echo "Refusing to start equilibrium production while on battery." >&2
  exit 75
fi

if pgrep -f '[e]ntropy_ft_affinity_equilibrium sample 6 6 10' >/dev/null; then
  echo "A matching equilibrium production is already active." >&2
  exit 76
fi

if ! mkdir "$HERE/.equilibrium.run.lock" 2>/dev/null; then
  echo "Equilibrium run lock already exists." >&2
  exit 77
fi

trap 'status=$?; printf "%s\n" "$status" > "$HERE/equilibrium_pipeline.exitcode"; exit "$status"' EXIT

if [[ -e "$BLOCKS" ]]; then
  echo "Refusing to overwrite existing equilibrium raw data: $BLOCKS" >&2
  exit 4
fi

actual_source_sha="$(shasum -a 256 "$SOURCE" | awk '{print $1}')"
if [[ "$actual_source_sha" != "$SOURCE_SHA" ]]; then
  echo "Frozen source hash mismatch: $actual_source_sha" >&2
  exit 2
fi

available_kb="$(df -Pk "$HERE" | awk 'NR==2 {print $4}')"
if [[ "$available_kb" -lt 1048576 ]]; then
  echo "At least 1 GiB free is required; found ${available_kb} KiB." >&2
  exit 5
fi

"$CXX" -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I"$EIGEN3_INCLUDE" \
  -I"$LIBOMP_PREFIX/include" \
  -L"$LIBOMP_PREFIX/lib" -lomp \
  "$SOURCE" -o "$BIN"

"$BIN" selftest

binary_sha="$(shasum -a 256 "$BIN" | awk '{print $1}')"
{
  date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
  printf 'experiment_repository_commit=%s\n' "$(git -C "$ROOT" rev-parse HEAD)"
  printf 'protocol_amendment=%s\n' 'AMENDMENT_001_EQUILIBRIUM.md'
  printf 'production_source_commit=%s\n' "$SOURCE_COMMIT"
  printf 'production_source_sha256=%s\n' "$SOURCE_SHA"
  printf 'binary_sha256=%s\n' "$binary_sha"
  printf 'parameters=T_left=6 T_right=6 n=10 gamma=0.1 dt=0.0005 burnin=500 block_time=20 batches=8 lanes=16 streams=128 blocks_per_stream=7813 blocks=1000064 threads=2 bond=5 seed=2026090405\n'
} > "$MANIFEST"

printf 'case\tcommand\n' > "$COMMANDS"
printf 'dbeta_0p000000\t%s sample 6 6 10 8 500 20 7813 0.0005 2026090405 2 %s 5\n' \
  "$BIN" "$OUTPUT_PREFIX" >> "$COMMANDS"

"$BIN" sample 6 6 10 8 500 20 7813 0.0005 2026090405 2 \
  "$OUTPUT_PREFIX" 5 > "$LOG" 2>&1 &
pid="$!"
printf '%s\n' "$pid" > "$HERE/equilibrium.pid"
wait "$pid"

expected_header='stream_id,block_id,q_left,q_right,delta_energy,entropy_medium,entropy_rate,action_current,energy_balance_error'
header="$(head -n 1 "$BLOCKS")"
if [[ "$header" != "$expected_header" ]]; then
  echo "Unexpected equilibrium CSV header." >&2
  exit 6
fi

rows="$(awk 'END {print NR-1}' "$BLOCKS")"
if [[ "$rows" != "1000064" ]]; then
  echo "Unexpected equilibrium row count: $rows" >&2
  exit 7
fi

if ! awk -F, 'NR > 1 {
  row = NR - 2;
  expected_stream = int(row / 7813);
  expected_block = row % 7813;
  if (NF != 9 || $1 + 0 != expected_stream || $2 + 0 != expected_block) exit 1;
  for (i=1; i<=NF; ++i) if ($i == "" || tolower($i) ~ /nan|inf/) exit 1;
} END {if (NR != 1000065) exit 1}' "$BLOCKS"; then
  echo "Equilibrium output finite/row/identifier audit failed." >&2
  exit 8
fi

hash="$(shasum -a 256 "$BLOCKS" | awk '{print $1}')"
{
  printf 'case\trows\tsha256\n'
  printf 'dbeta_0p000000\t%s\t%s\n' "$rows" "$hash"
} > "$RAW_MANIFEST"
date -u '+completed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$MANIFEST"
echo "Equilibrium n=10 production complete: $BLOCKS"
