# n=60 current pilot — 2026-06-20

## Material Passport

- Artifact type: code experiment result and validation note
- Model version: `gibbs-canonical-v1`
- Purpose: finite-size robustness check beyond the `n=50` larger-chain run.
- Status: completed production-size robustness extension; not used to redefine
  the manuscript's primary `n=10,20,30,40` exponent.

## Command

The burn-in continues the production/larger-chain rule
`burnin = 3.2 n^2`, giving `11520` for `n=60`.

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 60 8 11520 200 0.0005 20260623 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b8
```

Medium pilot:

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 60 16 11520 200 0.0005 20260624 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b16
```

Production-size robustness extension:

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 60 64 11520 200 0.0005 20260625 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b64
```

Matched production-resolution fine-timestep check:

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 60 64 11520 200 0.00025 20260626 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b64_dt2p5e-4
```

Diagnostic six-length fit including the production-size extension:

```sh
python3 flux/analyze_canonical_flux.py \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n10_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n20_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n30_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n40_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b64_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b64_dt2p5e-4_summary.csv \
  --primary-dt 0.0005 \
  --bootstrap 10000 \
  --seed 20260625 \
  --output-prefix Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n10_60_b64_scaling
```

## Result summary

Initial `n=60` pilot:

- trajectories: `128`
- burn-in: `11520`
- measurement window: `200`
- timestep: `5e-4`
- mean current: `0.01189574358`
- standard error: `0.00122139546`
- normal 95% CI: `[0.00950185247, 0.01428963468]`
- first-half/second-half stationarity statistic: `-0.447` paired SE
- elapsed time: `337.10` seconds

Medium `n=60` pilot:

- trajectories: `256`
- burn-in: `11520`
- measurement window: `200`
- timestep: `5e-4`
- mean current: `0.01361149053`
- standard error: `0.00081058642`
- normal 95% CI: `[0.01202277034, 0.01520021073]`
- first-half/second-half stationarity statistic: `0.350` paired SE
- elapsed time: `724.52` seconds

Production-size `n=60` run:

- trajectories: `1024`
- burn-in: `11520`
- measurement window: `200`
- timestep: `5e-4`
- mean current: `0.01244829643`
- standard error: `0.00041661977`
- normal 95% CI: `[0.01163173668, 0.01326485617]`
- first-half/second-half stationarity statistic: `-0.524` paired SE
- elapsed time: `9711.25` seconds

Matched production-resolution fine-timestep check:

- trajectories: `1024`
- burn-in: `11520`
- measurement window: `200`
- timestep: `2.5e-4`
- mean current: `0.01288458131`
- standard error: `0.00036914544`
- normal 95% CI: `[0.01216106954, 0.01360809308]`
- first-half/second-half stationarity statistic: `-0.866` paired SE
- elapsed time: `7298.92` seconds
- difference relative to the matched `dt=5e-4`, `1024`-trajectory run:
  `+0.00043628488`, or `+3.50%`, equal to `0.78` pooled standard errors

Six-length diagnostic fit using the four production lengths, the `n=50`
robustness run, and the `n=60` production-size extension:

- `E[J(n)] = 35.94 n^-1.930`
- log-fit `R^2 = 0.9974`
- bootstrap 95% exponent CI: `[-1.954,-1.906]`
- no stationarity flags with `|z| >= 2`

The production adjacent `n=50`--`60` local slope is `-2.178`.  The production
mean is about `10.5%` below the direct `n=10,20,30,40,50` extrapolation.  The
production mean differs from the medium pilot by `-1.28` pooled standard
errors and from the initial pilot by `0.43` pooled standard errors, so the
three `n=60` estimates are mutually compatible at the pilot/production Monte
Carlo level.

## Interpretation

The `n=60` production-size extension is consistent with the qualitative
faster-than-Fourier picture and does not indicate a crossover toward the
Fourier exponent over the available window.  The matched fine-timestep
production check is compatible with the `dt=5e-4` run within Monte Carlo
uncertainty.  Together these runs strengthen the manuscript's finite-size and
timestep robustness evidence but do not replace the primary `n=10,20,30,40`
exponent, because a systematic larger-length convergence study remains outside
the scope of the present revision.
