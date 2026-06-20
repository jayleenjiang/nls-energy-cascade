# Parameter-robustness pilot for action-current scaling

## Material Passport

- Artifact type: parameter-robustness pilot and production-upgrade plan
- Model version: `gibbs-canonical-v1`
- Origin date: 2026-06-20
- Verification status: PILOT ONLY
- Scope: bath-temperature robustness of the canonical action-current
  accumulator; not a manuscript-level numerical claim yet

## Motivation

The primary manuscript claim is intentionally scoped to the production
configuration `T1=10`, `Tn=2`, `gamma=0.1`, `dt=5e-4`, and chain lengths
`n=10,20,30,40`, with larger-chain robustness checks at `n=50,60`.

A natural reviewer question is whether the faster-than-Fourier finite-size
current decay is a peculiarity of this single bath-temperature choice.  This
pilot tests two bath-temperature perturbations while keeping the same canonical
source, timestep, burn-in rule, measurement window, and measured middle bond.
The source currently fixes `gamma=0.1`, so gamma-robustness is left as a
separate code-extension task.

## Pilot design

Common settings:

- simulator: `Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical`
- source: `flux/NLS_flux_canonical.cpp`
- `gamma=0.1`
- `dt=5e-4`
- measurement window: `200`
- chain lengths: `n=10,20,30,40`
- burn-ins: `1000,1280,2880,5120`
- pilot size: `4` batches = `64` trajectories per chain length
- threads: `2`

Parameter sets:

1. `moderate_contrast_T8_T4`: lower temperature contrast at comparable total
   bath scale.
2. `scaled_contrast_T5_T1`: same bath ratio as the primary run, but lower
   absolute bath scale.

This is deliberately a low-cost screening run.  The pilot is sufficient for
checking sign, rough scaling, stationarity red flags, and whether a
production-size upgrade is scientifically worthwhile.  It is not sufficient
for a final paper claim because the bootstrap intervals still reflect only
64 trajectories per length.

## Pilot commands

For each parameter set and chain length, the command form was

```sh
flux_canonical T1 Tn n 4 burnin 200 0.0005 seed 2 out_prefix
```

with seeds `20260626+n` for `T1=8,Tn=4` and `20260627+n` for `T1=5,Tn=1`.

The scaling fits were produced by

```sh
python3 flux/analyze_canonical_flux.py \
  .../n10_b4_summary.csv .../n20_b4_summary.csv \
  .../n30_b4_summary.csv .../n40_b4_summary.csv \
  --primary-dt 0.0005 \
  --bootstrap 10000 \
  --output-prefix .../b4_scaling
```

## Pilot results

| parameter set | fitted exponent | bootstrap 95% CI | log-fit `R^2` | max `|stationarity z|` |
|---|---:|---:|---:|---:|
| `T1=8,Tn=4` | `-1.8423` | `[-1.9724,-1.7281]` | `0.9899` | `1.18` |
| `T1=5,Tn=1` | `-1.8656` | `[-1.9742,-1.7689]` | `0.9979` | `1.97` |

Mean-current pilot table:

| parameter set | `n` | mean action current | SE | stationarity z |
|---|---:|---:|---:|---:|
| `T1=8,Tn=4` | 10 | `0.2277647330` | `0.0067477682` | `0.500` |
| `T1=8,Tn=4` | 20 | `0.0724420357` | `0.0031521335` | `-1.179` |
| `T1=8,Tn=4` | 30 | `0.0357905714` | `0.0023134629` | `-0.487` |
| `T1=8,Tn=4` | 40 | `0.0165741691` | `0.0016880753` | `0.447` |
| `T1=5,Tn=1` | 10 | `0.2006990119` | `0.0040343183` | `-1.974` |
| `T1=5,Tn=1` | 20 | `0.0612475350` | `0.0021057008` | `-0.475` |
| `T1=5,Tn=1` | 30 | `0.0258282586` | `0.0015487915` | `0.230` |
| `T1=5,Tn=1` | 40 | `0.0153969354` | `0.0013398923` | `-0.119` |

## Interpretation

Both pilot parameter sets show positive action currents and fitted exponents
well below `-1` over `n=10,20,30,40`.  This supports the decision to run a
production-size parameter-robustness check if the manuscript needs one more
reviewer-facing numerical reinforcement.

However, the result should remain outside the main manuscript claim until at
least one parameter set is upgraded to production scale.  The `T1=5,Tn=1`
pilot has a borderline first-half/second-half statistic at `n=10`
(`|z|=1.97`), just below the existing red-flag threshold `|z|>=2`; this is not
a failure, but it argues for production-size replication before using the
result in the paper.

## Recommended production upgrade

If the paper needs an additional robustness subsection, upgrade only one set
first:

1. `T1=8,Tn=4` is the cleaner first choice because the pilot stationarity
   diagnostics are comfortably below `|z|=2`.
2. Use the same production settings as the primary run:
   `64` batches = `1024` trajectories per length, `dt=5e-4`, measurement
   window `200`, and the existing `n`-dependent burn-ins.
3. After production completion, rerun the analyzer with `10000` bootstrap
   replicates.
4. Add a short manuscript paragraph only if all four production means are
   positive, no stationarity statistic has `|z|>=2`, and the fitted exponent
   remains below `-1` with a reasonable bootstrap interval.

Possible production commands:

```sh
flux_canonical 8 4 10 64 1000 200 0.0005 20260630 2 parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/n10
flux_canonical 8 4 20 64 1280 200 0.0005 20260630 2 parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/n20
flux_canonical 8 4 30 64 2880 200 0.0005 20260630 2 parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/n30
flux_canonical 8 4 40 64 5120 200 0.0005 20260630 2 parameter_robustness_2026-06-20/moderate_contrast_T8_T4_prod/n40
```

The production upgrade is optional.  The current manuscript is already locally
submission-gated without this extra parameter study; this artifact records the
next most useful numerical reinforcement if time and compute budget allow it.
