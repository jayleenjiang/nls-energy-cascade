#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
EXP="$ROOT/experiments/entropy_ft_n3_equilibrium_2026-09-01"

"$EXP/run_equilibrium_case.sh" 6 2026090106 T6
"$EXP/run_equilibrium_case.sh" 10 2026090110 T10
