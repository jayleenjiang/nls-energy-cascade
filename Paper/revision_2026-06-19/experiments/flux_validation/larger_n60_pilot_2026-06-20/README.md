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

Diagnostic six-length fit including this pilot:

```sh
python3 flux/analyze_canonical_flux.py \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n10_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n20_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n30_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n40_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n60_b8_summary.csv \
  --primary-dt 0.0005 \
  --bootstrap 10000 \
  --seed 20260623 \
  --output-prefix Paper/revision_2026-06-19/experiments/flux_validation/larger_n60_pilot_2026-06-20/n10_60_b8pilot_scaling
```

## Result summary

Pilot `n=60` run:

- trajectories: `128`
- burn-in: `11520`
- measurement window: `200`
- timestep: `5e-4`
- mean current: `0.01189574358`
- standard error: `0.00122139546`
- normal 95% CI: `[0.00950185247, 0.01428963468]`
- first-half/second-half stationarity statistic: `-0.447` paired SE
- elapsed time: `337.10` seconds

Six-length diagnostic fit using the four production lengths, the `n=50`
robustness run, and this `n=60` pilot:

- `E[J(n)] = 37.46 n^-1.944`
- log-fit `R^2 = 0.9965`
- bootstrap 95% exponent CI: `[-2.018,-1.885]`
- no stationarity flags with `|z| >= 2`

The adjacent `n=50`--`60` local slope is `-2.427`.  The pilot mean is about
`14.5%` below the direct `n=10,20,30,40,50` extrapolation, but this difference
is comparable to the pilot's wider Monte Carlo uncertainty.

## Interpretation

This pilot is consistent with the qualitative faster-than-Fourier picture and
does not indicate a crossover toward the Fourier exponent over the available
window.  It should not replace the manuscript's primary `n=10,20,30,40`
production exponent, because the `n=60` sample size is only `128` trajectories.
If the paper needs another numerical strengthening pass before submission, the
natural next step is a production-size `n=60` run, or a matched fine-timestep
pilot at the largest length.
