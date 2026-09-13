# Canonical flux production manifest

Run date: 2026-06-19  
Host: Apple M4, macOS 15.6.1  
Compiler: Apple clang 17.0.0  
Model version: `gibbs-canonical-v1`

## Frozen implementation

- Source: `flux/NLS_flux_canonical.cpp`
- Source SHA-256:
  `76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e`
- Binary:
  `Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical`
- Binary SHA-256:
  `1ebab166366fa475e3350cf60367c528669a52ea7380dc5b2572c5efc6da761d`

Build command:

```sh
clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/include/eigen3 \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  flux/NLS_flux_canonical.cpp \
  -o Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical
```

The build completed with no compiler warnings.

## Production configuration

Common parameters:

- `T1=10`, `Tn=2`, `gamma=0.1`
- fixed-step Euler--Maruyama, `dt=0.0005`
- 64 SIMD batches × 16 lanes = 1024 independent trajectories per length
- measurement window = 200
- base seed = 20260619
- two OpenMP threads per simultaneously launched run
- measured bond = `n/2`

Burn-in:

| n | burn-in |
|---:|---:|
| 10 | 1000 |
| 20 | 1280 |
| 30 | 2880 |
| 40 | 5120 |

Commands:

```sh
flux_canonical 10 2 10 64 1000 200 0.0005 20260619 2 production_dt5e-4/n10
flux_canonical 10 2 20 64 1280 200 0.0005 20260619 2 production_dt5e-4/n20
flux_canonical 10 2 30 64 2880 200 0.0005 20260619 2 production_dt5e-4/n30
flux_canonical 10 2 40 64 5120 200 0.0005 20260619 2 production_dt5e-4/n40
```

Each run records four equal current blocks, so the same trajectories support
finite-time averaging windows 50, 100, and 200.  The summary also reports a
paired first-half/second-half stationarity diagnostic and positivity-projection
frequency.

## Pre-production gates passed

- warning-clean optimized build;
- AddressSanitizer/UndefinedBehaviorSanitizer smoke run;
- byte-identical samples, profiles, and burn-in traces for one- versus
  two-thread runs with the same seed;
- equal-temperature mean current consistent with zero;
- swapping bath temperatures reverses the current with matching magnitude;
- pilot runs at all four chain lengths showed no `|stationarity z| >= 2`;
- `n=40` current at `dt=0.0005` and `dt=0.00025` differed by approximately
  0.03% in the 64-trajectory pilot.

## Production results

Primary analysis:

- `production_dt5e-4/flux_primary_runs.csv`
- `production_dt5e-4/flux_primary_scaling.json`
- `production_dt5e-4/flux_primary_scaling.pdf`

Current means:

| n | mean action current | SE | stationarity z |
|---:|---:|---:|---:|
| 10 | 0.3925219606 | 0.0018902080 | -0.4245 |
| 20 | 0.1191693526 | 0.0009195305 | 0.9021 |
| 30 | 0.0545731139 | 0.0006191849 | -0.7769 |
| 40 | 0.0297475540 | 0.0004827205 | -0.8856 |

Power-law fit:

- `E[J(n)] = 28.7457 n^-1.85008`
- log-fit `R^2 = 0.998013`
- trajectory-bootstrap exponent 95% CI: `[-1.87034, -1.83049]`
- production-only maximum first-half/second-half stationarity statistic:
  `|z| = 0.902`

Window analysis:

- `production_dt5e-4/current_windows_window_statistics.csv`
- `production_dt5e-4/current_windows_variance_scaling.pdf`
- per-`n` standardized finite-time-current histograms in
  `production_dt5e-4/current_windows_n*_windows.pdf`

Equal-temperature Gibbs validation:

- independent MCMC reference:
  `gibbs_mcmc/n6_T2_medium_profile.csv`
- SDE-vs-MCMC comparison:
  `gibbs_sde/gibbs_sde_n6_T2_summary.json`
- maximum SDE/MCMC profile discrepancy: `1.10` combined SE
- RMS relative profile difference: `0.269%`

Status: production run complete and analyzed.
