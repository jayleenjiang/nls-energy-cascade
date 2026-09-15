#!/usr/bin/env bash
set -euo pipefail

repo_root=$(cd "$(dirname "$0")/../.." && pwd)
source_file="$repo_root/flux/NLS_entropy_cloning.cpp"
binary="$repo_root/flux/entropy_cloning_v2"

clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  "$source_file" -o "$binary"

"$binary" selftest
