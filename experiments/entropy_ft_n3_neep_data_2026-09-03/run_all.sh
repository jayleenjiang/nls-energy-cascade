#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_neep_data_2026-09-03"

"$EXP/run_case.sh" 10 2 2026090310 driven
"$EXP/build_neep_archive.sh" driven
"$EXP/run_case.sh" 6 6 2026090306 equilibrium
"$EXP/build_neep_archive.sh" equilibrium
