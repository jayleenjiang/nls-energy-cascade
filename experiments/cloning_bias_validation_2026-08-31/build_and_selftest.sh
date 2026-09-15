#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
root="$repo_root/experiments/cloning_bias_validation_2026-08-31"
source_file="$repo_root/flux/NLS_entropy_cloning.cpp"
archive="$root/source_archive/NLS_entropy_cloning.cpp"
binary="$root/bin/entropy_cloning_v2"

mkdir -p "$root/source_archive" "$root/bin"
cp -p "$source_file" "$archive"

clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  "$archive" -o "$binary"

"$binary" selftest | tee "$root/SELFTEST.log"

{
  echo "# Build provenance"
  echo
  echo '```text'
  clang++ --version
  echo '```'
  echo
  echo "- Git HEAD: \`$(git -C "$repo_root" rev-parse HEAD)\`"
  echo "- source Git blob: \`$(git -C "$repo_root" hash-object "$source_file")\`"
  echo "- source SHA-256: \`$(shasum -a 256 "$archive" | awk '{print $1}')\`"
  echo "- executable SHA-256: \`$(shasum -a 256 "$binary" | awk '{print $1}')\`"
  echo "- build command: \`clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp source_archive/NLS_entropy_cloning.cpp -o bin/entropy_cloning_v2\`"
  echo "- self-test command: \`bin/entropy_cloning_v2 selftest\`"
} > "$root/BUILD_PROVENANCE.md"

shasum -a 256 \
  "$root/PROTOCOL.md" "$root/README.md" "$root/RUN_MATRIX.csv" \
  "$root/generate_run_matrix.py" "$root/build_and_selftest.sh" \
  "$root/run_validation.sh" "$root/run_pipeline_when_ac.sh" \
  "$root/analyze_bias_validation.py" "$archive" "$binary" \
  > "$root/PRE_RUN_MANIFEST.sha256"
