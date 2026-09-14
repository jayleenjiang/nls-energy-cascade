#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../../.." && pwd)"
RAW="$ROOT/experiments/spectral_gap_controlled_2026-09-07/raw"

find "$RAW/equilibrium_dt1e-3" "$RAW/equilibrium_dt2p5e-4" \
  -type f -name 'stream_*.f32' -exec stat -f '%b %N' {} + \
  | awk '$1 == 0 {$1=""; sub(/^ /, ""); print}' \
  | xargs -n 1 -P 8 sh -c \
      'dd if="$1" of=/dev/null bs=1m status=none && echo "HYDRATED $1"' sh

remaining="$(find "$RAW/equilibrium_dt1e-3" "$RAW/equilibrium_dt2p5e-4" \
  -type f -name 'stream_*.f32' -exec stat -f '%b %N' {} + \
  | awk '$1 == 0 {n++} END {print n+0}')"
echo "REMAINING_DATALESS $remaining"
test "$remaining" -eq 0
