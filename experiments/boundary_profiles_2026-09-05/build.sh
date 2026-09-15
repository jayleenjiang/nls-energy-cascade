#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
/opt/homebrew/bin/g++-15 -O3 -march=native -fopenmp -DNDEBUG -Wall -Wextra \
  -I/opt/homebrew/include/eigen3 "$ROOT/src/NLS_boundary_profiles.cpp" \
  -o "$ROOT/bin/NLS_boundary_profiles"
shasum -a 256 "$ROOT/source_archive/NLS_flux_SIMD_fixed.original.cpp" \
  "$ROOT/src/NLS_boundary_profiles.cpp" "$ROOT/bin/NLS_boundary_profiles" \
  > "$ROOT/BUILD_HASHES.sha256"
