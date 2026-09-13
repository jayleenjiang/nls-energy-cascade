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

Fine-timestep pilot:

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 50 16 8000 200 0.00025 20260622 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b16_dt2p5e-4
```

Production-resolution fine-timestep check:

```sh
Paper/revision_2026-06-19/experiments/flux_validation/bin/flux_canonical \
  10 2 50 64 8000 200 0.00025 20260623 4 \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_dt2p5e-4
```

Primary scaling-analysis command:

```sh
python3 flux/analyze_canonical_flux.py \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n10_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n20_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n30_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/production_dt5e-4/n40_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_summary.csv \
  Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n50_b64_dt2p5e-4_summary.csv \
  --primary-dt 0.0005 \
  --bootstrap 10000 \
  --seed 20260620 \
  --output-prefix Paper/revision_2026-06-19/experiments/flux_validation/larger_n_pilot_2026-06-20/n10_50_b64_scaling
```

Fit-window sensitivity command:

```sh
python3 Paper/revision_2026-06-19/scripts/analyze_flux_scaling_sensitivity.py
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

Fine-timestep pilot:

- configuration: `n=50`, `16` batches, `256` trajectories, burn-in `8000`,
  measurement window `200`, timestep `2.5e-4`
- mean current: `0.01879771710`
- standard error: `0.00081439495`
- normal 95% CI: `[0.01720153233, 0.02039390187]`
- first-half/second-half stationarity statistic: `-0.455` paired SE
- difference relative to the `dt=5e-4`, `1024`-trajectory run:
  `+0.00028187025`, or `+1.52%`, equal to `0.30` pooled standard errors

Production-resolution fine-timestep check:

- configuration: `n=50`, `64` batches, `1024` trajectories, burn-in `8000`,
  measurement window `200`, timestep `2.5e-4`
- mean current: `0.01918191598`
- standard error: `0.00040113161`
- normal 95% CI: `[0.01839571247, 0.01996811949]`
- first-half/second-half stationarity statistic: `-1.557` paired SE
- difference relative to the matched `dt=5e-4`, `1024`-trajectory run:
  `+0.00066606913`, or `+3.60%`, equal to `1.12` pooled standard errors

Fit-window sensitivity:

- primary `n=10,20,30,40`: exponent `-1.85008`, bootstrap 95% CI
  `[-1.87032,-1.83081]`, log-fit `R^2=0.99801`
- with `n=50`: exponent `-1.89449`, bootstrap 95% CI
  `[-1.91717,-1.87295]`, log-fit `R^2=0.99761`
- tail `n=20,30,40,50`: exponent `-2.03265`, bootstrap 95% CI
  `[-2.07868,-1.98781]`, log-fit `R^2=0.99935`
- adjacent local slopes range from `-1.71976` on `n=10--20` to
  `-2.12473` on `n=40--50`

## Interpretation

The larger-chain run is below the direct extrapolation of the original
`n=10,20,30,40` fit but remains consistent with the qualitative conclusion that
the current decays faster than the Fourier `1/n` scaling.  The
production-resolution fine-timestep check is consistent with the `dt=5e-4` run
within Monte Carlo error, with a `3.60%` shift equal to `1.12` pooled standard
errors.  The fit-window sensitivity analysis does not show drift toward the
Fourier exponent over the available window.  The manuscript should still keep
the `n=10,20,30,40` fit as the primary production exponent because `n=50` and
`n=60` are robustness extensions rather than a systematic larger-length
convergence study.
