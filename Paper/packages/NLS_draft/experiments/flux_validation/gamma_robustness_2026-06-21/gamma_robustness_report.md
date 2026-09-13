# Thermostat-coupling robustness report

Generated: `2026-06-21T02:33:22.393272+00:00`

Status: **PASS**

## Scope

This report tests whether the faster-than-Fourier finite-size decay of the
action current is tied to the single thermostat coupling used in the primary
production run.  Gamma-specific sources are generated from the frozen
canonical source; the frozen source itself is not edited.

## Frozen-source check

- Source: `flux/NLS_flux_canonical.cpp`
- Expected SHA-256: `76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e`
- Observed SHA-256: `76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e`
- Match: `True`

## Scaling results

| dataset | gamma | n values | exponent | 95% bootstrap CI | R^2 | max |z| |
|---|---:|---:|---:|---:|---:|---:|
| primary_reference_gamma0p1 | 0.1 | 10,20,30,40 | -1.85008 | [-1.87019, -1.83075] | 0.99801 | 0.90 |
| gamma0p05 | 0.05 | 10,20,30,40 | -1.65035 | [-1.66794, -1.63333] | 0.99382 | 1.14 |
| gamma0p2 | 0.2 | 10,20,30,40 | -1.99149 | [-2.01710, -1.96682] | 0.99963 | 1.74 |

## Per-length means

| dataset | gamma | n | mean action current | SE | stationarity z |
|---|---:|---:|---:|---:|---:|
| primary_reference_gamma0p1 | 0.1 | 10 | 0.3925219606 | 0.00189 | -0.42 |
| primary_reference_gamma0p1 | 0.1 | 20 | 0.1191693526 | 0.00092 | 0.90 |
| primary_reference_gamma0p1 | 0.1 | 30 | 0.05457311394 | 0.000619 | -0.78 |
| primary_reference_gamma0p1 | 0.1 | 40 | 0.02974755396 | 0.000483 | -0.89 |
| gamma0p05 | 0.05 | 10 | 0.3519908103 | 0.00166 | 1.14 |
| gamma0p05 | 0.05 | 20 | 0.1302638586 | 0.00093 | -1.00 |
| gamma0p05 | 0.05 | 30 | 0.06256684183 | 0.000641 | 0.62 |
| gamma0p05 | 0.05 | 40 | 0.03507030144 | 0.000494 | 0.55 |
| gamma0p2 | 0.2 | 10 | 0.3449211025 | 0.00177 | -1.74 |
| gamma0p2 | 0.2 | 20 | 0.090549237 | 0.000804 | -0.51 |
| gamma0p2 | 0.2 | 30 | 0.03971588599 | 0.000533 | -0.75 |
| gamma0p2 | 0.2 | 40 | 0.02167963165 | 0.00046 | 0.49 |

## Interpretation

All production-resolution thermostat-coupling datasets in this report are
interpreted only as finite-size robustness checks.  They support the narrower
claim that the observed faster-than-Fourier action-current decay over
`n=10,20,30,40` is not an artifact of the single `gamma=0.1` coupling.
They do not constitute an asymptotic transport theorem or a systematic
two-parameter bath sweep.

## Output files

- JSON: `Paper/revision_2026-06-19/experiments/flux_validation/gamma_robustness_2026-06-21/gamma_robustness_scaling.json`
- CSV: `Paper/revision_2026-06-19/experiments/flux_validation/gamma_robustness_2026-06-21/gamma_robustness_summary.csv`
- Plot PDF: `Paper/revision_2026-06-19/experiments/flux_validation/gamma_robustness_2026-06-21/gamma_robustness_scaling.pdf`
- Plot PNG: `Paper/revision_2026-06-19/experiments/flux_validation/gamma_robustness_2026-06-21/gamma_robustness_scaling.png`
