#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN_DIR="${1:-$ROOT/experiments/entropy_ft_2026-08-26/production}"
EXPECTED_BLOCKS="${2:-1000064}"
FINAL_DIR="${3:-$RUN_DIR/final_v1}"
STATUS_LABEL="${4:-production}"
ANALYSIS_DIR="${ANALYSIS_DIR:-$RUN_DIR/analysis}"
ADAPTIVE_DIR="$FINAL_DIR/adaptive"
VALIDATION_DIR="${VALIDATION_DIR:-$(dirname "$RUN_DIR")/validation}"
PYTHON="${PYTHON:-/opt/homebrew/bin/python3}"
LATEXMK="${LATEXMK:-}"
REQUIRE_HASH_MATCH="${REQUIRE_HASH_MATCH:-1}"

if [[ -e "$FINAL_DIR" ]]; then
  echo "Refusing to overwrite existing finalization directory: $FINAL_DIR" >&2
  exit 2
fi
if [[ ! -d "$ANALYSIS_DIR" ]]; then
  echo "Core analysis directory does not exist: $ANALYSIS_DIR" >&2
  exit 3
fi
manifest="$(find "$RUN_DIR" -maxdepth 1 -name '*manifest.txt' -print -quit)"
if [[ -z "$manifest" ]] || ! grep -q '^completed_utc=' "$manifest"; then
  echo "Run manifest is absent or not marked complete: ${manifest:-none}" >&2
  exit 4
fi
if grep -q '^failed_utc=' "$manifest"; then
  echo "Run manifest contains failed_utc; refusing finalization." >&2
  exit 5
fi

mkdir -p "$FINAL_DIR/supplement" "$ADAPTIVE_DIR" "$FINAL_DIR/report"
{
  date -u '+started_utc=%Y-%m-%dT%H:%M:%SZ'
  printf 'run_dir=%s\n' "$RUN_DIR"
  printf 'analysis_dir=%s\n' "$ANALYSIS_DIR"
  printf 'validation_dir=%s\n' "$VALIDATION_DIR"
  printf 'expected_blocks_per_n=%s\n' "$EXPECTED_BLOCKS"
  printf 'require_hash_match=%s\n' "$REQUIRE_HASH_MATCH"
  git -C "$ROOT" rev-parse HEAD
  shasum -a 256 \
    "$ROOT/flux/audit_entropy_ft_run.py" \
    "$ROOT/flux/plot_entropy_ft_supplement.py" \
    "$ROOT/flux/analyze_ft_adaptive_bins.py" \
    "$ROOT/flux/audit_entropy_ft_analysis.py" \
    "$ROOT/flux/build_entropy_ft_report.py"
} > "$FINAL_DIR/finalization_manifest.txt"

raw_audit_args=(
  "$ROOT/flux/audit_entropy_ft_run.py"
  "$RUN_DIR"
  --expected-blocks "$EXPECTED_BLOCKS"
  --require-completed
  --output-prefix "$FINAL_DIR/raw_audit"
)
if [[ "$REQUIRE_HASH_MATCH" == 1 ]]; then
  raw_audit_args+=(--require-hash-match)
fi
"$PYTHON" "${raw_audit_args[@]}"

"$PYTHON" "$ROOT/flux/plot_entropy_ft_supplement.py" \
  "$RUN_DIR"/n*_blocks.csv \
  --analysis-dir "$ANALYSIS_DIR" \
  --output-dir "$FINAL_DIR/supplement" \
  --taus 20,40,60,80,100,120,140,160,180,200 \
  --threshold-step 0.01 \
  --minimum-raw-count 5 \
  --tail-probability-max 0.01

"$PYTHON" "$ROOT/flux/analyze_ft_adaptive_bins.py" \
  "$RUN_DIR"/n*_blocks.csv \
  --output-dir "$ADAPTIVE_DIR" \
  --taus 20,40,60,80,100,120,140,160,180,200 \
  --max-bins 60 \
  --min-effective-count 50 \
  --range-quantile 0.99 \
  --bootstrap 1000

"$PYTHON" "$ROOT/flux/audit_entropy_ft_analysis.py" \
  "$RUN_DIR" \
  --analysis-dir "$ANALYSIS_DIR" \
  --supplement-dir "$FINAL_DIR/supplement" \
  --adaptive-dir "$ADAPTIVE_DIR" \
  --output-prefix "$FINAL_DIR/analysis_audit"

"$PYTHON" "$ROOT/flux/build_entropy_ft_report.py" \
  --run-dir "$RUN_DIR" \
  --analysis-dir "$ANALYSIS_DIR" \
  --supplement-dir "$FINAL_DIR/supplement" \
  --adaptive-dir "$ADAPTIVE_DIR" \
  --validation-dir "$VALIDATION_DIR" \
  --output-dir "$FINAL_DIR/report" \
  --status-label "$STATUS_LABEL"

if [[ -z "$LATEXMK" ]]; then
  if command -v latexmk >/dev/null 2>&1; then
    LATEXMK="$(command -v latexmk)"
  elif [[ -x /Library/TeX/texbin/latexmk ]]; then
    LATEXMK=/Library/TeX/texbin/latexmk
  else
    echo "latexmk is required to compile the final report." >&2
    exit 6
  fi
fi
mkdir -p "$FINAL_DIR/report/build"
(
  cd "$FINAL_DIR/report"
  "$LATEXMK" -norc -pdf -interaction=nonstopmode -halt-on-error \
    -outdir=build entropy_ft_report.tex
)

if grep -Eq 'Overfull|Underfull|undefined references|LaTeX Warning' \
  "$FINAL_DIR/report/build/entropy_ft_report.log"; then
  echo "Final report log contains a layout or reference warning." >&2
  exit 7
fi

shasum -a 256 "$RUN_DIR"/n*_blocks.csv > "$FINAL_DIR/raw_blocks_sha256.txt"
date -u '+completed_utc=%Y-%m-%dT%H:%M:%SZ' \
  >> "$FINAL_DIR/finalization_manifest.txt"
{
  find "$RUN_DIR" -maxdepth 1 -type f ! -name 'n*_blocks.csv' -print0
  find "$ANALYSIS_DIR" -type f -print0
  find "$FINAL_DIR" -type f \
    ! -name 'curated_sha256.txt' \
    ! -path "$FINAL_DIR/report/build/*" \
    -print0
} | sort -zu \
  | xargs -0 shasum -a 256 > "$FINAL_DIR/curated_sha256.txt"
shasum -a 256 "$FINAL_DIR/report/build/entropy_ft_report.pdf" \
  >> "$FINAL_DIR/curated_sha256.txt"

echo "Finalization complete: $FINAL_DIR"
