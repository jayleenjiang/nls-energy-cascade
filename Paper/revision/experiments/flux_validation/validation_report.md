# Flux experiment validation report

## Material Passport

- Origin Skill: experiment-agent
- Origin Mode: validate
- Origin Date: 2026-06-19
- Verification Status: VERIFIED
- Version Label: canonical_flux_validation_v1

## Scope

This report validates the corrected long-chain action-current experiment for
the Gibbs-preserving NLS SDE. It replaces the June 18 flux data, which used
boundary noise amplitudes missing the required factor `sqrt(2)`.

## Reproducibility and implementation checks

- Frozen source SHA-256:
  `76f937608280272397a555931b353ba770b06ee87d2f5b0dce08fe1e6bb3727e`
- Frozen binary SHA-256:
  `1ebab166366fa475e3350cf60367c528669a52ea7380dc5b2572c5efc6da761d`
- Build: Apple clang 17.0.0, warning-clean.
- Sanitizer smoke run: AddressSanitizer/UndefinedBehaviorSanitizer completed.
- Determinism: same seed produced byte-identical samples, profiles, and
  burn-in traces with one and two OpenMP threads.
- Physical smoke tests:
  - Equal-temperature run had zero action current within CI.
  - Swapping hot/cold baths reversed the current with matching magnitude.
  - Equal-temperature SDE action profile agreed with an independent
    Metropolis sample of `exp(-H/(2T))`; maximum discrepancy was `1.10`
    combined SE.

## Primary statistical findings

| n | mean action current | SE | 95% CI | stationarity z |
|---:|---:|---:|---:|---:|
| 10 | 0.3925219606 | 0.0018902080 | [0.3888172210, 0.3962267003] | -0.4245 |
| 20 | 0.1191693526 | 0.0009195305 | [0.1173671060, 0.1209715992] | 0.9021 |
| 30 | 0.0545731139 | 0.0006191849 | [0.0533595338, 0.0557866940] | -0.7769 |
| 40 | 0.0297475540 | 0.0004827205 | [0.0288014392, 0.0306936687] | -0.8856 |

Scaling fit over `n = 10,20,30,40`:

- `E[J(n)] = 28.7457 n^-1.85008`
- log-fit `R^2 = 0.998013`
- trajectory-bootstrap exponent 95% CI: `[-1.87034, -1.83049]`

Larger-chain robustness and fit-window sensitivity:

- the `n=50` robustness run gives `E[J(50)] = 0.01851584685` with SE
  `0.00044158954`
- the `n=60` production-size robustness run gives
  `E[J(60)] = 0.01244829643` with SE `0.00041661977`
- adding `n=50` gives a diagnostic exponent `-1.89449`, bootstrap 95% CI
  `[-1.91717, -1.87295]`, and log-fit `R^2 = 0.99761`
- adding `n=50,60` gives a diagnostic exponent `-1.92956`, bootstrap 95%
  CI `[-1.95424, -1.90603]`, and log-fit `R^2 = 0.99739`
- fitting only the tail `n=20,30,40,50,60` gives exponent `-2.05926`,
  bootstrap 95% CI `[-2.10697, -2.01242]`, and log-fit `R^2 = 0.99933`
- adjacent local slopes range from `-1.71976` on `n=10--20` to `-2.17771`
  on `n=50--60`

Bath-temperature robustness:

- a second production-resolution four-length run at `T1=8,Tn=4`, `gamma=0.1`,
  and `dt=5e-4` gives positive currents at all four primary lengths
- the means are `0.2233162853`, `0.07169996295`, `0.03422040656`, and
  `0.01948996976` for `n=10,20,30,40`
- the corresponding standard errors are `0.00177600849`, `0.00086742549`,
  `0.00056584124`, and `0.00044790010`
- the fitted finite-size exponent is `-1.75098`, bootstrap 95% CI
  `[-1.77964, -1.72269]`, with log-fit `R^2 = 0.99844`
- the maximum first-half/second-half stationarity statistic is `1.73684`
  paired standard errors, below the `|z|>=2` red-flag threshold
- this is a bath-temperature robustness check, not a systematic parameter
  sweep

Timestep sensitivity:

- `dt=1e-3` is visibly coarse at larger `n` and underestimates the
  `n=40` current.
- `dt=5e-4` and `dt=2.5e-4` agree within Monte Carlo error in the pilot:
  for `n=40`, relative difference `3.9%`, `0.60` pooled SE.

Finite-window current distributions:

- Windows `tau=50,100,200` are produced from four disjoint current blocks per
  trajectory.
- At `n=40`, `Pr(Jbar_tau<0)` is `0.233`, `0.118`, and `0.026` for the three
  windows.
- The product `tau Var(Jbar_tau)` decreases from `0.095` to `0.048`, so the
  data are treated as finite-window descriptive statistics, not as an
  asymptotic large-deviation estimate.

## Statistical fallacy scan

Coverage: 11/11 checked.

| Fallacy | Status | Note |
|---|---|---|
| Simpson's paradox | Not applicable | No grouped aggregate reversal analysis. |
| Ecological fallacy | Not applicable | Inference stays at trajectory/chain-level simulation. |
| Berkson's paradox | Not detected | No selected conditional sample. |
| Collider bias | Not detected | No regression with controls. |
| Base-rate neglect | Not applicable | No diagnostic/screening probabilities. |
| Regression to the mean | Not applicable | No selected extreme pre/post groups. |
| Survivorship bias | Not detected | All trajectories are retained; no dropouts. |
| Look-elsewhere effect | Caution addressed | Exponent chosen after correcting model; CI reported; no selective p-values. |
| Garden of forking paths | Caution addressed | Pilot/production split and timestep decision recorded in manifest. |
| Correlation != causation | Not applicable | Mechanistic SDE simulation, not observational causal claim. |
| Reverse causality | Not applicable | No causal inference from cross-sectional data. |

## Reproducibility verdict

VERIFIED for the canonical action-current experiment. The evidence supports the
manuscript claim that, for the corrected Gibbs-preserving SDE at
`T1=10,Tn=2,gamma=0.1`, the finite-chain mean action current decays
approximately as `n^-1.85` over `n=10,20,30,40`. A production-resolution
robustness check at `T1=8,Tn=4,gamma=0.1` gives a comparable faster-than-Fourier
finite-size decay with exponent approximately `-1.75`.

Limitations: the exponent is a finite-size numerical scaling over four chain
lengths; the confidence interval reflects Monte Carlo uncertainty conditional
on this scaling model and does not by itself prove an asymptotic theorem.
