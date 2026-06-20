# n=60 current pilot — 2026-06-20

## Material Passport

- Artifact type: code experiment result and pilot validation note
- Model version: `gibbs-canonical-v1`
- Purpose: feasibility and finite-size robustness check beyond the `n=50`
  larger-chain run.
- Status: completed pilot; not used as a primary manuscript claim.

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

Diagnostic six-length fit including the medium pilot:

```sh
python3 flux/analyze_canonical_flux.py \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n10_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n20_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n30_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n40_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b16_summary.csv \
  --primary-dt 0.0005 \
  --bootstrap 10000 \
  --seed 20260624 \
  --output-prefix Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n10_60_b16pilot_scaling
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

Six-length diagnostic fit using the four production lengths, the `n=50`
robustness run, and the `n=60` medium pilot:

- `E[J(n)] = 33.14 n^-1.901`
- log-fit `R^2 = 0.9982`
- bootstrap 95% exponent CI: `[-1.943,-1.864]`
- no stationarity flags with `|z| >= 2`

The medium-pilot adjacent `n=50`--`60` local slope is `-1.688`.  The medium
pilot mean is about `2.1%` below the direct `n=10,20,30,40,50` extrapolation.
The initial and medium pilots differ by `1.17` pooled standard errors, so their
spread is compatible with pilot-scale Monte Carlo variation.

## Interpretation

Both `n=60` pilots are consistent with the qualitative faster-than-Fourier
picture and do not indicate a crossover toward the Fourier exponent over the
available window.  They should not replace the manuscript's primary
`n=10,20,30,40` production exponent, because the larger `n=60` sample size is
still only `256` trajectories.  If the paper needs another numerical
strengthening pass before submission, the natural next step is a production-size
`n=60` run, or a matched fine-timestep pilot at the largest length.
