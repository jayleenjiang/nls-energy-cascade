# Reproduction guide

Run commands from the repository root:

```sh
cd /Users/jayleenjiang/Documents/NLS
```

The production environment was Apple Silicon, Clang with OpenMP, and Python
3.  Exact Phase-II environment notes are in
`results/long_chain_phase2/ENVIRONMENT.md`.

## 1. Direct two-tail production

Build and self-test the projection-free Cartesian sampler:

```sh
clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  flux/NLS_entropy_ft.cpp -o flux/entropy_ft

./flux/entropy_ft selftest
```

The exact production parameters and source hashes are in
`protocols/direct_sampling_production_manifest.txt`.  The original launch
wrapper is `code/direct_sampling/run_entropy_ft_production.sh`.

Re-run the fail-fast finalization and all audits:

```sh
flux/finalize_entropy_ft_run.sh
```

The wrapper verifies exactly 1,000,064 rows for each chain length, recomputes
heat, entropy and first-law columns, performs the fixed-bin and adaptive-bin
analyses, builds both-tail survival curves, fits the descriptive normal-tail
benchmark, analyses `log P` versus averaging time, and writes the report.

## 2. Two-site total entropy

The accepted endpoint-density analysis is reproduced by:

```sh
python3 flux/analyze_total_entropy_parametric_n2.py \
  experiments/entropy_ft_scgf_2026-08-27/total_entropy_n2_short/production/equilibrium_blocks.csv \
  experiments/entropy_ft_scgf_2026-08-27/total_entropy_n2_short/production/driven_blocks.csv \
  --output-dir /tmp/n2_parametric_reproduction

python3 flux/audit_total_entropy_parametric_n2.py \
  /tmp/n2_parametric_reproduction \
  --output /tmp/n2_parametric_reproduction/audit.md
```

The frozen acceptance conditions are in
`results/n2_total_entropy/PREDECLARED_GATES.md`.

## 3. Exact discrete path-ratio control

```sh
clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  flux/NLS_discrete_path_ft.cpp -o flux/discrete_path_ft

./flux/discrete_path_ft selftest

python3 flux/analyze_discrete_path_ft.py \
  experiments/discrete_path_ft_2026-08-28/production_v2/driven_t0p1_dt1e3_N1m_forward.csv \
  experiments/discrete_path_ft_2026-08-28/production_v2/driven_t0p1_dt1e3_N1m_reverse.csv \
  --output-dir /tmp/discrete_path_reproduction
```

Only the `v2` reverse implementation is accepted.  Earlier exploratory
reverse files used the Hamiltonian sign twice and are retained only for
provenance.

## 4. Long-chain SCGF Phase II

```sh
experiments/long_chain_ft_phase2_2026-08-29/build_and_selftest.sh
experiments/long_chain_ft_phase2_2026-08-29/analyze_phase2.sh final
```

The simulation matrix and all gates were frozen in `PROTOCOL.md`.  The package
contains all 56 Phase-II summaries, timeseries, and logs, so final analysis can
be reproduced without the large direct-sampling CSV files.

## 5. Compile the combined report

```sh
python3 /Users/jayleenjiang/.codex/plugins/cache/openai-bundled/latex/0.2.6/scripts/compile_latex.py \
  /Users/jayleenjiang/Documents/NLS/Paper/packages/mentor_round2_ft_2026-08-30/mentor_round2_ft_report.tex \
  --compiler texlive \
  --output-directory /Users/jayleenjiang/Documents/NLS/Paper/packages/mentor_round2_ft_2026-08-30/build
```
