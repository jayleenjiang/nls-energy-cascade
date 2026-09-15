#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
BASE="$HERE/../entropy_ft_window_length_2026-09-10"
BIN="$BASE/bin/entropy_ft_window"
SOURCE="$BASE/source/NLS_entropy_ft_1905cf.cpp"
RAW="$HERE/raw"
PYTHON="${PYTHON:-/opt/homebrew/bin/python3}"
SOURCE_SHA="98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3"
BINARY_SHA="4c4880d721733897d200f2601690da873f5b226aaa1013df23df2351c1dfe7d1"
EXPECTED_HEADER='stream_id,block_id,q_left,q_right,delta_energy,entropy_medium,entropy_rate,action_current,energy_balance_error'
EXPECTED_ROWS=4000256
BLOCKS_PER_STREAM=31252

mkdir -p "$RAW" "$HERE/analysis" "$HERE/figures" "$HERE/report"
if ! pmset -g batt | head -1 | grep -q 'AC Power'; then
  echo 'Refusing to start production while on battery.' >&2
  exit 75
fi
if pgrep -f '[e]ntropy_ft_window sample' >/dev/null; then
  echo 'A matching entropy production is already active.' >&2
  exit 76
fi
if ! mkdir "$HERE/.run.lock" 2>/dev/null; then
  echo "Run lock already exists: $HERE/.run.lock" >&2
  exit 77
fi
trap 'status=$?; printf "%s\n" "$status" > "$HERE/pipeline.exitcode"; exit "$status"' EXIT

[[ "$(shasum -a 256 "$SOURCE" | awk '{print $1}')" == "$SOURCE_SHA" ]]
[[ "$(shasum -a 256 "$BIN" | awk '{print $1}')" == "$BINARY_SHA" ]]
"$BIN" selftest

available_kb="$(df -Pk "$HERE" | awk 'NR==2 {print $4}')"
if [[ "$available_kb" -lt 8388608 ]]; then
  echo "At least 8 GiB free is required; found ${available_kb} KiB." >&2
  exit 78
fi

cases=(weak_dt2p5e4 weak_dt1p25e4 moderate_dt2p5e4 moderate_dt1p25e4)
tls=(6.5 6.5 8 8)
trs=(5.5 5.5 4 4)
dts=(0.00025 0.000125 0.00025 0.000125)
seeds=(2026091101 2026091102 2026091103 2026091104)

for case_name in "${cases[@]}"; do
  output="$RAW/${case_name}_n10_tau5_blocks.csv"
  if [[ -e "$output" ]]; then
    echo "Refusing to overwrite $output" >&2
    exit 79
  fi
done

{
  date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
  printf 'repository_commit_at_launch=%s\n' "$(git -C "$ROOT" rev-parse HEAD)"
  printf 'production_source_commit=%s\n' '1905cf4e606a4a7f4dd8930caa64bd4cc861e9d4'
  printf 'source_path=%s\n' "$SOURCE"
  printf 'source_sha256=%s\n' "$SOURCE_SHA"
  printf 'binary_path=%s\n' "$BIN"
  printf 'binary_sha256=%s\n' "$BINARY_SHA"
  printf 'parameters=n=10 gamma=0.1 burnin=500 base_tau=5 blocks_per_stream=31252 streams=128 rows_per_case=4000256 stationary_time_per_stream=156260 total_stationary_time_per_case=20001280 threads_per_case=2 bond=5\n'
} > "$HERE/production_manifest.txt"

printf 'case\tcommand\n' > "$HERE/COMMANDS.tsv"
printf 'case\tpid\tstart_epoch\n' > "$HERE/production_pids.tsv"
: > "$HERE/power_events.log"

pids=()
starts=()
for i in 0 1 2 3; do
  case_name="${cases[$i]}"
  prefix="$RAW/${case_name}_n10_tau5"
  printf '%s\t%s sample %s %s 10 8 500 5 31252 %s %s 2 %s 5\n' \
    "$case_name" "$BIN" "${tls[$i]}" "${trs[$i]}" "${dts[$i]}" \
    "${seeds[$i]}" "$prefix" >> "$HERE/COMMANDS.tsv"
  start="$(date +%s)"
  (exec "$BIN" sample "${tls[$i]}" "${trs[$i]}" 10 8 500 5 31252 \
    "${dts[$i]}" "${seeds[$i]}" 2 "$prefix" 5) \
    > "$RAW/${case_name}_n10_tau5.log" 2>&1 &
  pids+=("$!")
  starts+=("$start")
  printf '%s\t%s\t%s\n' "$case_name" "$!" "$start" >> "$HERE/production_pids.tsv"
done

done_flags=(0 0 0 0)
printf 'case\tstart_epoch\tend_epoch\telapsed_seconds\n' > "$HERE/RUNTIMES.tsv"
stopped=0
remaining=4
while [[ "$remaining" -gt 0 ]]; do
  for i in 0 1 2 3; do
    [[ "${done_flags[$i]}" -eq 1 ]] && continue
    pid="${pids[$i]}"
    if ! kill -0 "$pid" 2>/dev/null; then
      end="$(date +%s)"
      printf '%s\t%s\t%s\t%s\n' "${cases[$i]}" "${starts[$i]}" "$end" \
        "$((end - starts[$i]))" >> "$HERE/RUNTIMES.tsv"
      done_flags[$i]=1
      remaining="$((remaining - 1))"
    fi
  done
  [[ "$remaining" -eq 0 ]] && break
  if pmset -g batt | head -1 | grep -q 'AC Power'; then
    if [[ "$stopped" -eq 1 ]]; then
      for i in 0 1 2 3; do
        [[ "${done_flags[$i]}" -eq 0 ]] && kill -CONT "${pids[$i]}" 2>/dev/null || true
      done
      date -u '+%Y-%m-%dT%H:%M:%SZ resumed_on_ac' >> "$HERE/power_events.log"
      stopped=0
    fi
  elif [[ "$stopped" -eq 0 ]]; then
    for i in 0 1 2 3; do
      [[ "${done_flags[$i]}" -eq 0 ]] && kill -STOP "${pids[$i]}" 2>/dev/null || true
    done
    date -u '+%Y-%m-%dT%H:%M:%SZ stopped_on_battery' >> "$HERE/power_events.log"
    stopped=1
  fi
  sleep 30
done

status=0
for pid in "${pids[@]}"; do
  if ! wait "$pid"; then status=1; fi
done
if [[ "$status" -ne 0 ]]; then
  echo 'At least one simulator failed; preserving all outputs.' >&2
  exit "$status"
fi

printf 'case\trows\tsha256\n' > "$HERE/RAW_DATA_MANIFEST.tsv"
for case_name in "${cases[@]}"; do
  blocks="$RAW/${case_name}_n10_tau5_blocks.csv"
  [[ "$(head -n 1 "$blocks")" == "$EXPECTED_HEADER" ]]
  rows="$(awk 'END {print NR-1}' "$blocks")"
  [[ "$rows" == "$EXPECTED_ROWS" ]]
  awk -F, -v bps="$BLOCKS_PER_STREAM" 'NR>1 {
    row=NR-2; stream=int(row/bps); block=row%bps;
    if (NF!=9 || $1+0!=stream || $2+0!=block) exit 1;
    for (i=1;i<=NF;i++) if ($i=="" || tolower($i) ~ /nan|inf/) exit 1;
  } END {if (NR-1 != 4000256) exit 1}' "$blocks"
  hash="$(shasum -a 256 "$blocks" | awk '{print $1}')"
  printf '%s\t%s\t%s\n' "$case_name" "$rows" "$hash" >> "$HERE/RAW_DATA_MANIFEST.tsv"
done

"$PYTHON" "$HERE/analyze_multiaffinity_timestep.py" --experiment-dir "$HERE"
"$PYTHON" "$HERE/build_report.py" --experiment-dir "$HERE"
(
  cd "$HERE/report"
  /Library/TeX/texbin/latexmk -pdf -interaction=nonstopmode -halt-on-error \
    multiaffinity_timestep_report.tex
)

date -u '+completed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$HERE/production_manifest.txt"
echo "Multi-affinity timestep production and analysis complete: $HERE"

