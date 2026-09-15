# Fixed SIMD flux distribution rerun

Run timestamp: 2026-06-19 22:40 EDT

## Purpose

Minimal debug/rerun of the original SIMD flux-distribution workflow.  The
original June 18 `NLS_flux_SIMD.cpp` output format is preserved: each run writes
`flux_n<N>.txt` containing one header line followed by one finite-time averaged
flux sample per trajectory, plus an appended `flux_vs_length_2.csv` summary.

## Source files

- Original source: `flux/flux_data/NLS_flux_SIMD.cpp`
- Fixed source: `flux/flux_data/NLS_flux_SIMD_fixed.cpp`
- Fixed binary: `flux/flux_data/NLS_flux_SIMD_fixed`

SHA-256:

```text
be188aef7d644f23be153db614be73a52a7d31b798f809a2870f473494ec863b  flux/flux_data/NLS_flux_SIMD.cpp
14ab6f2cf70255c458d29e6d4686b81b212e8af150e93a9fb4a945f909b84bd4  flux/flux_data/NLS_flux_SIMD_fixed.cpp
9e1dc7260357203d8d856852b0700925d6f991fc364a038542d1d3893fc74373  flux/flux_data/NLS_flux_SIMD_fixed
```

## Minimal fixes

1. Boundary noise amplitudes now match the Gibbs-preserving SDE:
   `2*sqrt(2*gamma*T*I)` for action noise and `sqrt(2*gamma*T/I)` for phase
   noise.
2. RNG seeding is deterministic by default with base seed `20260619`; seed
   depends only on the SIMD batch index, not on OpenMP thread assignment.
3. Printed and CSV confidence intervals use `1.96 * std_dev / sqrt(N)` rather
   than `1.96 * std_dev`.

The implementation otherwise keeps the original SIMD/float/adaptive-step
workflow and output names.

## Build command

```sh
clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/include/eigen3 \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  flux/flux_data/NLS_flux_SIMD_fixed.cpp \
  -o flux/flux_data/NLS_flux_SIMD_fixed
```

The fixed source compiled with no warnings.

## Reproducibility smoke test

Two independent 2-batch `n=10` smoke runs with the same seed produced
byte-identical `flux_n10.txt` files (`cmp` exit code `0`).  The smoke runs were
executed under `/private/tmp/nls_simd_repro_a` and
`/private/tmp/nls_simd_repro_b` and did not touch the production outputs in this
directory.

## Run commands

All runs use `case=0`, `T1=10`, `Tn=2`, `batches=625` (`10000` trajectories),
`T_burnin=1000`, `T_final=1200`, measurement window `200`, and base seed
`20260619`.

```sh
../NLS_flux_SIMD_fixed 0 10 2 10 625 1200 20260619 > n10.log 2>&1
../NLS_flux_SIMD_fixed 0 10 2 20 625 1200 20260619 > n20.log 2>&1
../NLS_flux_SIMD_fixed 0 10 2 30 625 1200 20260619 > n30.log 2>&1
../NLS_flux_SIMD_fixed 0 10 2 40 625 1200 20260619 > n40.log 2>&1
```

## Outputs

- Raw samples: `flux_n10.txt`, `flux_n20.txt`, `flux_n30.txt`, `flux_n40.txt`
- SIMD summary: `flux_vs_length_2.csv`
- Tail analysis: `simd_fixed_tail_summary.csv`, `flux_distribution.png`
- Canonical comparison: `simd_fixed_vs_canonical.csv`,
  `simd_fixed_scaling_compare.json`, `simd_fixed_vs_canonical.png`,
  `simd_fixed_vs_canonical.pdf`
