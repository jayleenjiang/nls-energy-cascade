#!/usr/bin/env bash
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../.." && pwd)"
OUT="${1:-$HOME/Downloads/nls_multiaffinity_timestep_audit_2026-09-13.tar.zst}"

export COPYFILE_DISABLE=1
tar \
  --exclude='*.aux' \
  --exclude='*.fdb_latexmk' \
  --exclude='*.fls' \
  --exclude='*.out' \
  --exclude='__pycache__' \
  --exclude='*.pyc' \
  --exclude='*.pid' \
  --exclude='.run.lock' \
  -C "$REPO" -cf - \
  experiments/entropy_ft_multiaffinity_timestep_2026-09-11 \
  experiments/entropy_ft_window_length_2026-09-10/raw/dbeta_0p027972_n10_tau5_blocks.csv \
  experiments/entropy_ft_window_length_2026-09-10/raw/dbeta_0p125000_n10_tau5_blocks.csv \
  experiments/entropy_ft_window_length_2026-09-10/raw/dbeta_0p000000_n10_tau5_blocks.csv \
  experiments/entropy_ft_window_length_2026-09-10/source/NLS_entropy_ft_1905cf.cpp \
  experiments/entropy_ft_window_length_2026-09-10/bin/entropy_ft_window \
  experiments/entropy_ft_window_length_2026-09-10/analyze_window_dependence.py \
  experiments/entropy_ft_timestep_bias_2026-09-10/raw/equilibrium_dt2p5e4_n10_tau5_blocks.csv \
  experiments/entropy_ft_timestep_bias_2026-09-10/analysis/extrapolation_summary.csv \
  experiments/entropy_ft_pooled_extrapolation_2026-09-10/analyze_pooled_extrapolation.py \
  | zstd -T0 -6 -f -o "$OUT"

zstd -t "$OUT"
shasum -a 256 "$OUT" > "$OUT.sha256"
printf 'package=%s\n' "$OUT"
printf 'sha256_file=%s.sha256\n' "$OUT"
