#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
SOURCE="$HERE/source/NLS_entropy_ft_1905cf.cpp"
BIN="$HERE/bin/entropy_ft_affinity"
RAW="$HERE/raw"
ANALYSIS="$HERE/analysis"
FIGURES="$HERE/figures"
REPORT="$HERE/report"
MANIFEST="$HERE/production_manifest.txt"
COMMANDS="$HERE/COMMANDS.tsv"
PYTHON="${PYTHON:-/opt/homebrew/bin/python3}"
CXX="${CXX:-clang++}"
EIGEN3_INCLUDE="${EIGEN3_INCLUDE:-/opt/homebrew/include/eigen3}"
LIBOMP_PREFIX="${LIBOMP_PREFIX:-/opt/homebrew/opt/libomp}"
EXISTING_DRIVEN="/Users/jayleenjiang/Downloads/n10_blocks.csv"
SOURCE_COMMIT="1905cf4e606a4a7f4dd8930caa64bd4cc861e9d4"
SOURCE_SHA="98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3"
DRIVEN_SHA="a23806e82f5514a9c3375d10a6644946b6b57d0efe8452b3bbf397b6230f9929"

mkdir -p "$HERE/bin" "$RAW" "$ANALYSIS" "$FIGURES" "$REPORT"

if ! pmset -g batt | head -1 | grep -q "AC Power"; then
  echo "Refusing to start production while the Mac is on battery power." >&2
  exit 75
fi

if ! mkdir "$HERE/.run.lock" 2>/dev/null; then
  echo "Run lock already exists: $HERE/.run.lock" >&2
  exit 76
fi
trap 'status=$?; printf "%s\n" "$status" > "$HERE/pipeline.exitcode"; exit "$status"' EXIT

actual_source_sha="$(shasum -a 256 "$SOURCE" | awk '{print $1}')"
if [[ "$actual_source_sha" != "$SOURCE_SHA" ]]; then
  echo "Frozen source hash mismatch: $actual_source_sha" >&2
  exit 2
fi

actual_driven_sha="$(shasum -a 256 "$EXISTING_DRIVEN" | awk '{print $1}')"
if [[ "$actual_driven_sha" != "$DRIVEN_SHA" ]]; then
  echo "Existing (10,2) input hash mismatch: $actual_driven_sha" >&2
  exit 3
fi

cases=(dbeta_0p027972 dbeta_0p057143 dbeta_0p125000 dbeta_0p222222)
tls=(6.5 7 8 9)
trs=(5.5 5 4 3)
seeds=(2026090401 2026090402 2026090403 2026090404)

for case_name in "${cases[@]}"; do
  if [[ -e "$RAW/${case_name}_n10_blocks.csv" ]]; then
    echo "Refusing to overwrite existing raw data: $RAW/${case_name}_n10_blocks.csv" >&2
    exit 4
  fi
done

available_kb="$(df -Pk "$HERE" | awk 'NR==2 {print $4}')"
if [[ "$available_kb" -lt 2097152 ]]; then
  echo "At least 2 GiB free is required; found ${available_kb} KiB." >&2
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
analysis_sha="$(shasum -a 256 "$HERE/analyze_affinity_sweep.py" | awk '{print $1}')"
report_builder_sha="$(shasum -a 256 "$HERE/build_report.py" | awk '{print $1}')"

{
  date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
  printf 'experiment_repository_commit=%s\n' "$(git -C "$ROOT" rev-parse HEAD)"
  printf 'production_source_commit=%s\n' "$SOURCE_COMMIT"
  printf 'production_source_sha256=%s\n' "$SOURCE_SHA"
  printf 'binary_sha256=%s\n' "$binary_sha"
  printf 'analysis_sha256=%s\n' "$analysis_sha"
  printf 'report_builder_sha256=%s\n' "$report_builder_sha"
  printf 'existing_TL10_TR2_sha256=%s\n' "$DRIVEN_SHA"
  printf 'parameters=n=10 gamma=0.1 dt=0.0005 burnin=500 block_time=20 batches=8 lanes=16 streams=128 blocks_per_stream=7813 blocks_per_case=1000064 threads_per_case=2 bond=5\n'
  printf 'equilibrium_n10_status=UNAVAILABLE_NOT_RERUN_PER_PROTOCOL\n'
} > "$MANIFEST"

printf 'case\tcommand\n' > "$COMMANDS"
pids=()
for i in 0 1 2 3; do
  case_name="${cases[$i]}"
  tl="${tls[$i]}"
  tr="${trs[$i]}"
  seed="${seeds[$i]}"
  prefix="$RAW/${case_name}_n10"
  printf '%s\t%s sample %s %s 10 8 500 20 7813 0.0005 %s 2 %s 5\n' \
    "$case_name" "$BIN" "$tl" "$tr" "$seed" "$prefix" >> "$COMMANDS"
  "$BIN" sample "$tl" "$tr" 10 8 500 20 7813 0.0005 "$seed" 2 \
    "$prefix" 5 > "$RAW/${case_name}_n10.log" 2>&1 &
  pids+=("$!")
  printf '%s\t%s\n' "$case_name" "$!" >> "$HERE/production_pids.tsv"
done

status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then
    status=1
  fi
done
if [[ "$status" -ne 0 ]]; then
  date -u '+failed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$MANIFEST"
  echo "At least one affinity case failed; raw outputs are preserved." >&2
  exit "$status"
fi

printf 'case\trows\tsha256\n' > "$HERE/RAW_DATA_MANIFEST.tsv"
for case_name in "${cases[@]}"; do
  blocks="$RAW/${case_name}_n10_blocks.csv"
  header="$(head -n 1 "$blocks")"
  expected_header='stream_id,block_id,q_left,q_right,delta_energy,entropy_medium,entropy_rate,action_current,energy_balance_error'
  if [[ "$header" != "$expected_header" ]]; then
    echo "Unexpected header in $blocks" >&2
    exit 6
  fi
  rows="$(awk 'END {print NR-1}' "$blocks")"
  if [[ "$rows" != "1000064" ]]; then
    echo "Unexpected row count in $blocks: $rows" >&2
    exit 7
  fi
  hash="$(shasum -a 256 "$blocks" | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "$case_name" "$rows" "$hash" >> "$HERE/RAW_DATA_MANIFEST.tsv"
done

"$PYTHON" "$HERE/analyze_affinity_sweep.py" \
  --case "dbeta_0p027972:6.5:5.5:$RAW/dbeta_0p027972_n10_blocks.csv" \
  --case "dbeta_0p057143:7:5:$RAW/dbeta_0p057143_n10_blocks.csv" \
  --case "dbeta_0p125000:8:4:$RAW/dbeta_0p125000_n10_blocks.csv" \
  --case "dbeta_0p222222:9:3:$RAW/dbeta_0p222222_n10_blocks.csv" \
  --case "dbeta_0p400000:10:2:$EXISTING_DRIVEN" \
  --output-dir "$ANALYSIS" \
  --figures-dir "$FIGURES" \
  --bootstrap 1000 \
  --seed 2026090499

"$PYTHON" "$HERE/build_report.py" \
  --experiment-dir "$HERE" \
  --output "$REPORT/affinity_sweep_report.tex"

(
  cd "$REPORT"
  /Library/TeX/texbin/latexmk -pdf -interaction=nonstopmode -halt-on-error \
    affinity_sweep_report.tex
)

date -u '+completed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$MANIFEST"
echo "Affinity sweep production and analysis complete: $HERE"
