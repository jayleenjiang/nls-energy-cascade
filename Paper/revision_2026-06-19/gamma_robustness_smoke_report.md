# Gamma robustness smoke report

Generated: `2026-06-20T21:31:49.465402+00:00`

Status: **PASS**

## Scope

This is a code-path smoke test for future thermostat-coupling robustness runs.
It is not manuscript evidence and is not used in the current action-current scaling claim.

The frozen production source is not edited.  The script verifies its SHA-256,
generates gamma-specific temporary sources under `tmp/`, compiles them, and
runs tiny `n=6` simulations only to check the build/run/output path.

## Frozen-source check

- Source: `flux/NLS_flux_canonical.cpp`
- Expected SHA-256: `76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e`
- Observed SHA-256: `76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e`
- Match: `True`

## Smoke results

| gamma | compile | run | summary gamma | mean current | SE | notes |
|---:|---:|---:|---:|---:|---:|---|
| 0.05 | 0 | 0 | 0.05 | 0.5326551657790191 | 0.2892532714180685 | smoke ok |
| 0.2 | 0 | 0 | 0.2 | 0.611765378601208 | 0.5580905327281221 | smoke ok |

## Next production step

If the authors want to include gamma robustness as manuscript evidence,
run production-resolution chains for one or two gamma values, for example
`gamma=0.05` and `gamma=0.2`, using the same primary lengths
`n=10,20,30,40`, the validated timestep, the existing burn-in schedule,
and trajectory-level bootstrap analysis.  Do not fold pilot results into
the primary exponent unless the full production protocol passes the same
stationarity and metadata gates as the existing current-scaling data.
