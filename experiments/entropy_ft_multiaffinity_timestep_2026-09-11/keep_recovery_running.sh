#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PID_FILE="$HERE/recovery_pids.tsv"
LOG="$HERE/unconditional_continuation.log"

[[ -r "$PID_FILE" ]]
date -u '+%Y-%m-%dT%H:%M:%SZ unconditional-continuation monitor started' >> "$LOG"

while true; do
  alive=0
  while IFS=$'\t' read -r case_name pid start_epoch; do
    [[ "$case_name" == "case" ]] && continue
    [[ "$pid" =~ ^[0-9]+$ ]] || continue
    if kill -0 "$pid" 2>/dev/null; then
      alive=$((alive + 1))
      state="$(ps -o state= -p "$pid" | tr -d '[:space:]')"
      if [[ "$state" == *T* ]]; then
        kill -CONT "$pid"
        date -u "+%Y-%m-%dT%H:%M:%SZ SIGCONT case=$case_name pid=$pid" >> "$LOG"
      fi
    fi
  done < "$PID_FILE"
  [[ "$alive" -eq 0 ]] && break
  sleep 5
done

date -u '+%Y-%m-%dT%H:%M:%SZ unconditional-continuation monitor finished' >> "$LOG"
