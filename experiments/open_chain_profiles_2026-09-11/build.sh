#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "$0")" && pwd)"
mkdir -p "$ROOT/bin"
/opt/homebrew/bin/g++-15 -O3 -march=native -fopenmp -DNDEBUG -Wall -Wextra \
  -I/opt/homebrew/include/eigen3 "$ROOT/src/NLS_open_chain_profiles.cpp" \
  -o "$ROOT/bin/NLS_open_chain_profiles"
shasum -a 256 "$ROOT/src/NLS_open_chain_profiles.cpp" \
  "$ROOT/bin/NLS_open_chain_profiles" > "$ROOT/BUILD_HASHES.sha256"
