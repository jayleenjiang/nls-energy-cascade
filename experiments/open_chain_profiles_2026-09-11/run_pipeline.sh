#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BIN="$HERE/bin/NLS_open_chain_profiles"
MATRIX="$HERE/RUN_MATRIX.csv"
PYTHON="${PYTHON:-/opt/homebrew/bin/python3}"
SOURCE_SHA="dfee7ed0b416a30f2b27d919c524dfaf37e91ca44312ffbb9d68b7fc1cd36110"
BINARY_SHA="09d801e71b7f3bd13cea56aba3d432b471dbf16a726a9866894148a91397575d"

mkdir -p "$HERE/raw" "$HERE/logs" "$HERE/analysis" "$HERE/figures" \
  "$HERE/curated_results" "$HERE/report" "$HERE/package"
if pgrep -f '[N]LS_open_chain_profiles 4' >/dev/null; then
  echo 'Matching open-chain profile process already active; refusing duplicate.' >&2
  exit 76
fi
if ! mkdir "$HERE/.run.lock" 2>/dev/null; then
  echo "Run lock already exists: $HERE/.run.lock" >&2
  exit 77
fi
trap 'status=$?; printf "%s\n" "$status" > "$HERE/pipeline.exitcode"; exit "$status"' EXIT

[[ "$(shasum -a 256 "$HERE/src/NLS_open_chain_profiles.cpp" | awk '{print $1}')" == "$SOURCE_SHA" ]]
[[ "$(shasum -a 256 "$BIN" | awk '{print $1}')" == "$BINARY_SHA" ]]
printf 'case,command\n' > "$HERE/COMMANDS.csv"
printf 'case,start_epoch,end_epoch,elapsed_seconds\n' > "$HERE/RUNTIMES.csv"
: > "$HERE/power_events.log"
{
  date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
  printf 'repository_commit_at_launch=%s\n' "$(git -C "$HERE" rev-parse HEAD)"
  printf 'remote_freeze_commit=%s\n' "$(tr -d '[:space:]' < "$HERE/FREEZE_COMMIT.txt")"
  printf 'source_sha256=%s\n' "$SOURCE_SHA"
  printf 'binary_sha256=%s\n' "$BINARY_SHA"
  printf 'parent_profile_source_sha256=%s\n' '40d13c7547d9c147ca8241ed701e1a1e85540fc0c21db48c4a63a943eccfe3ec'
} > "$HERE/production_manifest.txt"

wait_for_ac() {
  while ! pmset -g batt | head -1 | grep -q 'AC Power'; do
    sleep 30
  done
}

run_one() {
  local label="$1" tl="$2" tr="$3" n="$4" rep="$5" batches="$6"
  local burnin="$7" measure="$8" seed="$9" threads="${10}"
  local prefix="$HERE/raw/OPEN_${label}_n${n}_rep${rep}"
  local log="$HERE/logs/OPEN_${label}_n${n}_rep${rep}.log"
  local required=(profile quarter_profiles quarter_differences burnin_checkpoints \
    checkpoints state_checkpoints mass_checkpoints trajectory_diagnostics summary)
  local complete=1 suffix
  for suffix in "${required[@]}"; do
    [[ -s "${prefix}_${suffix}.csv" ]] || complete=0
  done
  if [[ "$complete" -eq 1 ]]; then
    echo "skip_complete $prefix"
    return
  fi
  for suffix in "${required[@]}"; do
    if [[ -e "${prefix}_${suffix}.csv" ]]; then
      echo "Partial output exists; refusing overwrite: ${prefix}_${suffix}.csv" >&2
      exit 79
    fi
  done
  wait_for_ac
  local case_name="OPEN_${label}_n${n}_rep${rep}"
  local command="$BIN 4 $tl $tr $n $batches $burnin $measure $seed $prefix $threads"
  printf '%s,"%s"\n' "$case_name" "$command" >> "$HERE/COMMANDS.csv"
  local start end stopped=0
  start="$(date +%s)"
  "$BIN" 4 "$tl" "$tr" "$n" "$batches" "$burnin" "$measure" \
    "$seed" "$prefix" "$threads" > "$log" 2>&1 &
  local pid=$!
  printf '%s\n' "$pid" > "$HERE/active_simulator.pid"
  while kill -0 "$pid" 2>/dev/null; do
    if pmset -g batt | head -1 | grep -q 'AC Power'; then
      if [[ "$stopped" -eq 1 ]]; then
        kill -CONT "$pid" 2>/dev/null || true
        date -u "+%Y-%m-%dT%H:%M:%SZ resumed pid=$pid case=$case_name" >> "$HERE/power_events.log"
        stopped=0
      fi
    elif [[ "$stopped" -eq 0 ]]; then
      kill -STOP "$pid" 2>/dev/null || true
      date -u "+%Y-%m-%dT%H:%M:%SZ stopped pid=$pid case=$case_name" >> "$HERE/power_events.log"
      stopped=1
    fi
    sleep 30
  done
  if ! wait "$pid"; then
    echo "Simulator failed: $case_name" >&2
    exit 80
  fi
  end="$(date +%s)"
  printf '%s,%s,%s,%s\n' "$case_name" "$start" "$end" "$((end-start))" >> "$HERE/RUNTIMES.csv"
  for suffix in "${required[@]}"; do test -s "${prefix}_${suffix}.csv"; done
}

while IFS=, read -r label tl tr n rep batches burnin measure seed threads; do
  [[ "$label" == "temperature_label" ]] && continue
  run_one "$label" "$tl" "$tr" "$n" "$rep" "$batches" "$burnin" \
    "$measure" "$seed" "$threads"
done < "$MATRIX"

printf 'path,sha256\n' > "$HERE/RAW_HASHES.csv"
find "$HERE/raw" -type f -name '*.csv' | LC_ALL=C sort | while IFS= read -r path; do
  printf '"%s",%s\n' "$path" "$(shasum -a 256 "$path" | awk '{print $1}')" >> "$HERE/RAW_HASHES.csv"
done

"$PYTHON" "$HERE/analyze_open_chain.py" --experiment-dir "$HERE"
"$PYTHON" "$HERE/build_report.py" --experiment-dir "$HERE"
(
  cd "$HERE/report"
  /Library/TeX/texbin/latexmk -pdf -interaction=nonstopmode -halt-on-error \
    open_chain_profiles_report.tex
)

date -u '+completed_utc=%Y-%m-%dT%H:%M:%SZ' >> "$HERE/production_manifest.txt"
(
  cd "$HERE"
  /usr/bin/zip -q -r package/open_chain_profiles_code_results.zip \
    PROTOCOL.md RUN_MATRIX.csv src bin build.sh analyze_open_chain.py \
    build_report.py run_pipeline.sh run_when_ac.sh COMMANDS.csv RUNTIMES.csv \
    FREEZE_COMMIT.txt production_manifest.txt RAW_HASHES.csv raw analysis curated_results figures \
    report/open_chain_profiles_report.tex report/open_chain_profiles_report.pdf
)

echo 'Open-chain profile production and analysis complete.'
