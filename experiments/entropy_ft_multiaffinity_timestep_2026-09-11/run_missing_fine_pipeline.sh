#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
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
  echo 'Refusing to start recovery while on battery.' >&2
  exit 75
fi
if pgrep -f '[e]ntropy_ft_window sample' >/dev/null; then
  echo 'A matching entropy production is already active.' >&2
  exit 76
fi
if ! mkdir "$HERE/.fine_recovery.lock" 2>/dev/null; then
  echo "Recovery lock already exists: $HERE/.fine_recovery.lock" >&2
  exit 77
fi
trap 'status=$?; printf "%s\n" "$status" > "$HERE/recovery.exitcode"; exit "$status"' EXIT

[[ "$(shasum -a 256 "$SOURCE" | awk '{print $1}')" == "$SOURCE_SHA" ]]
[[ "$(shasum -a 256 "$BIN" | awk '{print $1}')" == "$BINARY_SHA" ]]
"$BIN" selftest

for case_name in weak_dt2p5e4 moderate_dt2p5e4; do
  blocks="$RAW/${case_name}_n10_tau5_blocks.csv"
  [[ -s "$blocks" ]]
done
for case_name in weak_dt1p25e4 moderate_dt1p25e4; do
  blocks="$RAW/${case_name}_n10_tau5_blocks.csv"
  if [[ -e "$blocks" ]]; then
    echo "Refusing to overwrite $blocks" >&2
    exit 79
  fi
done

available_kb="$(df -Pk "$HERE" | awk 'NR==2 {print $4}')"
if [[ "$available_kb" -lt 4194304 ]]; then
  echo "At least 4 GiB free is required; found ${available_kb} KiB." >&2
  exit 78
fi

{
  date -u '+recovery_started_utc=%Y-%m-%dT%H:%M:%SZ'
  printf 'reason=original fine-timestep processes exited before output; cause unresolved\n'
  printf 'source_sha256=%s\n' "$SOURCE_SHA"
  printf 'binary_sha256=%s\n' "$BINARY_SHA"
  printf 'parameters_unchanged=n=10 gamma=0.1 burnin=500 base_tau=5 blocks_per_stream=31252 streams=128 rows_per_case=4000256 threads_per_case=2 bond=5\n'
} > "$HERE/recovery_manifest.txt"

cases=(weak_dt1p25e4 moderate_dt1p25e4)
tls=(6.5 8)
trs=(5.5 4)
seeds=(2026091102 2026091104)
pids=()
starts=()
printf 'case\tcommand\n' > "$HERE/RECOVERY_COMMANDS.tsv"
printf 'case\tpid\tstart_epoch\n' > "$HERE/recovery_pids.tsv"
: > "$HERE/recovery_power_events.log"

for i in 0 1; do
  case_name="${cases[$i]}"
  prefix="$RAW/${case_name}_n10_tau5"
  printf '%s\t%s sample %s %s 10 8 500 5 31252 0.000125 %s 2 %s 5\n' \
    "$case_name" "$BIN" "${tls[$i]}" "${trs[$i]}" "${seeds[$i]}" "$prefix" \
    >> "$HERE/RECOVERY_COMMANDS.tsv"
  start="$(date +%s)"
  (exec "$BIN" sample "${tls[$i]}" "${trs[$i]}" 10 8 500 5 31252 \
    0.000125 "${seeds[$i]}" 2 "$prefix" 5) \
    > "$RAW/${case_name}_n10_tau5.recovery.log" 2>&1 &
  pids+=("$!")
  starts+=("$start")
  printf '%s\t%s\t%s\n' "$case_name" "$!" "$start" >> "$HERE/recovery_pids.tsv"
done

done_flags=(0 0)
statuses=(0 0)
printf 'case\tstart_epoch\tend_epoch\telapsed_seconds\texit_status\n' > "$HERE/RECOVERY_RUNTIMES.tsv"
stopped=0
remaining=2
while [[ "$remaining" -gt 0 ]]; do
  for i in 0 1; do
    [[ "${done_flags[$i]}" -eq 1 ]] && continue
    pid="${pids[$i]}"
    if ! kill -0 "$pid" 2>/dev/null; then
      if wait "$pid"; then statuses[$i]=0; else statuses[$i]=$?; fi
      end="$(date +%s)"
      printf '%s\t%s\t%s\t%s\t%s\n' "${cases[$i]}" "${starts[$i]}" "$end" \
        "$((end - starts[$i]))" "${statuses[$i]}" >> "$HERE/RECOVERY_RUNTIMES.tsv"
      done_flags[$i]=1
      remaining="$((remaining - 1))"
    fi
  done
  [[ "$remaining" -eq 0 ]] && break
  if pmset -g batt | head -1 | grep -q 'AC Power'; then
    if [[ "$stopped" -eq 1 ]]; then
      for i in 0 1; do
        [[ "${done_flags[$i]}" -eq 0 ]] && kill -CONT "${pids[$i]}" 2>/dev/null || true
      done
      date -u '+%Y-%m-%dT%H:%M:%SZ resumed_on_ac' >> "$HERE/recovery_power_events.log"
      stopped=0
    fi
  elif [[ "$stopped" -eq 0 ]]; then
    for i in 0 1; do
      [[ "${done_flags[$i]}" -eq 0 ]] && kill -STOP "${pids[$i]}" 2>/dev/null || true
    done
    date -u '+%Y-%m-%dT%H:%M:%SZ stopped_on_battery' >> "$HERE/recovery_power_events.log"
    stopped=1
  fi
  sleep 30
done

for status in "${statuses[@]}"; do
  if [[ "$status" -ne 0 ]]; then
    echo 'At least one recovery simulator failed; preserving all outputs.' >&2
    exit 1
  fi
done

all_cases=(weak_dt2p5e4 weak_dt1p25e4 moderate_dt2p5e4 moderate_dt1p25e4)
printf 'case\trows\tsha256\n' > "$HERE/RAW_DATA_MANIFEST.tsv"
for case_name in "${all_cases[@]}"; do
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

date -u '+recovery_completed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$HERE/recovery_manifest.txt"
echo "Fine-timestep recovery, audit, and analysis complete: $HERE"
