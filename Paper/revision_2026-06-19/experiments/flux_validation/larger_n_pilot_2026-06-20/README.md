# Larger-chain current robustness check — 2026-06-20

## Material Passport

- Artifact type: code experiment result and validation note
- Model version: `gibbs-canonical-v1`
- Purpose: test whether the corrected action-current scaling over
  `n=10,20,30,40` remains compatible with a larger chain length.
- Status: completed; suitable as a robustness check, not a replacement for a
  full larger-length/timestep convergence study.

## Commands

Primary larger-chain run:

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 50 64 8000 200 0.0005 20260620 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64
```

Initial smoke/pilot run:

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 50 8 8000 200 0.0005 20260620 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b8
```

Longer-burn-in sanity check:

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 50 16 10000 200 0.0005 20260621 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b16_burn10000
```

Primary scaling-analysis command:

```sh
python3 flux/analyze_canonical_flux.py \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n10_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n20_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n30_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n40_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_summary.csv \
  --primary-dt 0.0005 \
  --bootstrap 10000 \
  --seed 20260620 \
  --output-prefix Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n10_50_b64_scaling
```

## Result summary

Primary `n=50` run:

- trajectories: `1024`
- burn-in: `8000`
- measurement window: `200`
- timestep: `5e-4`
- mean current: `0.01851584685`
- standard error: `0.00044158954`
- normal 95% CI: `[0.01765034726, 0.01938134644]`
- first-half/second-half stationarity statistic: `-1.485` paired SE

Five-length diagnostic fit using the four production lengths plus this
larger-chain run:

- `E[J(n)] = 32.50 n^-1.894`
- log-fit `R^2 = 0.9976`
- bootstrap 95% exponent CI: `[-1.916, -1.873]`

Longer-burn-in sanity check:

- configuration: `n=50`, `16` batches, `256` trajectories, burn-in `10000`
- mean current: `0.01931242054`
- standard error: `0.00085798987`
- normal 95% CI: `[0.01763079128, 0.02099404979]`

## Interpretation

The larger-chain run is below the direct extrapolation of the original
`n=10,20,30,40` fit but remains consistent with the qualitative conclusion that
the current decays faster than the Fourier `1/n` scaling.  Because the `n=50`
run has not yet been paired with a fine-timestep sensitivity check, the
manuscript should keep the `n=10,20,30,40` fit as the primary production
exponent and use `n=50` as a robustness check.
