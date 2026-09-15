# Final verdict: multi-affinity timestep audit

Final status: **MODERATE_AFFINITY_ASYMPTOTE_UNRESOLVED**.

All four new simulations completed and passed the integrity audit.  Each raw
block file contains exactly 4,000,256 finite rows from 128 independent streams
with contiguous stream and block identifiers.  The frozen analysis was rerun
from the raw inputs and reproduced every numerical CSV and PNG exactly.

## Predeclared primary gate at the finest timestep

| Affinity | Delta beta | dt | a_inf / Delta beta | Bootstrap 95% CI | 1/tau fit adequate? | CI contains 1? | Verdict |
|---|---:|---:|---:|---:|---:|---:|---|
| (6.5,5.5) | 0.0279720 | 1.25e-4 | 0.993519 | [0.983582, 1.001209] | yes | yes | consistent with FT |
| (7,5), existing control | 0.0571429 | 1.25e-4 | 0.995730 | [0.988535, 1.000882] | yes | yes | consistent with FT |
| (8,4) | 0.125000 | 1.25e-4 | 0.959884 | [0.952289, 0.965475] | **no** | no | **asymptote unresolved** |

For `(8,4)`, the primary fit uses `tau=5,10,20,25,40,80` and fails both
frozen adequacy diagnostics: reduced chi-square `5.10276 > 2` and maximum
absolute standardized residual `3.24993 > 3`.  Therefore its primary
intercept is not accepted as an asymptotic FT estimate.  The independently
predeclared `tau>=10` fit is adequate at the finest step, but gives
`a_inf/Delta beta=0.978279` with 95% CI `[0.966215,0.987892]`, which excludes
one.  This strengthens the need for caution; it does not override the primary
gate.

## Timestep and estimator-floor interpretation

- Weak affinity: the diagnostic `dt -> 0` ratio is `0.995263`, 95% CI
  `[0.983447,1.004505]`, and is interpretable under the frozen gates.
- Moderate affinity: the diagnostic `dt -> 0` ratios are not interpretable
  because the constituent finite-time primary models are not all adequate.
- The equilibrium matched-bin estimator floor is `1.08828e-4` in absolute
  slope units.  Both equilibrium baselines contain zero, and their difference
  CI contains zero.
- The mean-heat timestep shifts are affinity-dependent and non-monotone; no
  universal correction sign or magnitude is supported.

## Claim boundary

The data support a paper-level statement that the heat-flux fluctuation
relation is numerically consistent at the weak `(6.5,5.5)` and intermediate
`(7,5)` affinities at the finest tested timestep.  They do **not** support a
three-affinity verification claim.  At `(8,4)`, the accessible negative-tail
support and the frozen finite-time model do not yield a reliable asymptotic
intercept; the result must be reported as unresolved, not as verified.
