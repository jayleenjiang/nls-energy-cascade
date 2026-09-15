# Predeclared gates for long-chain GC verification

These gates are fixed before the multi-seed, multi-population production
comparison is inspected.  Pilot runs may select the proposal control scale and
selection interval, but production acceptance uses the frozen values.

## Observable and estimator identity

1. Heat is accumulated separately at the two baths and is positive into the
   system.
2. The primary observable is
   `Sigma_R = -(1/T_R-1/T_L) Q_R`.
3. The recorded identity
   `Sigma_R-Sigma_medium=(1/T_L) Delta E + residual`
   must hold to the first-law numerical tolerance.
4. At zero control, the exact-controlled implementation must reproduce the
   naive path-weight implementation to machine precision for identical seeds.
5. Every production run must have zero nonfinite trajectories and zero
   midpoint failures.

## Independent-run and support gates

For every accepted `(n,t,k,N_c,dt,selection)` group:

- at least four independent seeds;
- minimum selection-weight ESS at least `0.1 N_c` in every run;
- finite SCGF and late-time log-normalizer slope in every run.

At the `t=20` calibration horizon, additionally require at least 32 surviving
initial roots and minimum root-weight ESS at least 16 in every run.  Initial
root diversity is reported at all later times but is not a hard long-time gate:
at fixed population it necessarily decays under repeated resampling even for a
correct Feynman--Kac model.  Long-time bias is instead controlled by independent
runs, per-selection ESS, population doubling, and late-time slope convergence.
This distinction was fixed after the adaptive-resampling pilot and before the
multi-seed long-time production runs.

## Paired fluctuation-symmetry gate

Only nontrivial pairs `k` and `1-k` are evidence; `k=1/2` is excluded.
For an accepted pair:

- both members separately pass every support gate;
- the independent-seed confidence interval for
  `D_n,t(k)=psi_n,t(k)-psi_n,t(1-k)` contains zero;
- `|mean D_n,t(k)| <= 0.01`;
- the conclusion is unchanged when using late-half log-normalizer slopes.

At least the pair `(0.3,0.7)` must pass at every chain length.  The pair
`(0.4,0.6)` is included whenever both endpoints pass support.

## Convergence gates

1. **Population:** doubling `N_c` changes each paired SCGF by no more than two
   combined independent-run SE and `0.01` in absolute value.
2. **Observation time:** the accepted residual is checked at no fewer than
   three cumulative horizons.  The final two horizons must differ by no more
   than two combined SE and `0.01`; they must retain full support.
3. **Timestep:** halving `dt` changes each member and its paired residual by no
   more than two combined SE and `0.01`.
4. **Selection interval:** at least two intervals must agree within two
   combined SE and `0.01`, with full support.
5. **Gauge:** the long-time medium-entropy and right-current SCGFs must approach
   one another where both direct or controlled estimates have support.  A
   finite-time difference is reported as an endpoint effect, not hidden.

## Scope of an accepted claim

Passing these gates supports numerical consistency with the long-time
Gallavotti--Cohen symmetry over the resolved `k` interval and tested chain
lengths.  It is not a mathematical proof of ergodicity or a statement for all
chain lengths and all tilts.

## Post-pilot amendment and claim scope

The original strongest acceptance target above required `(0.3,0.7)` at every
chain length.  The production support audit found that this pair is resolved at
`n=10` but that its `k=0.7` member is rare-event/finite-time limited for
`n>=20`.  Those failed diagnostics are retained.  They do **not** satisfy the
original all-chain strongest-pair gate.

The pair `(0.4,0.6)` was already predeclared as a secondary pair.  It is used to
state a narrower result for `n=20,30,40` only where all remaining gates pass.
Accordingly, the admissible conclusion is numerical consistency with GC over
the *resolved* pair at each tested chain length, together with an explicit
statement that `(0.3,0.7)` is unresolved on the longer chains.  It is not a
claim that the original strongest all-chain target passed.

The support-only selection-interval audit also replaced
`selection_time=0.5` by `selection_time=2.0` before the replacement production
series.  The earlier interval and all associated output remain pilot evidence;
the reason for the amendment is documented in `PILOT_DECISIONS.md`.

## Post-production n=40 remedial extension

The four-seed `n=40`, `N_c=1024` series passed support and full-window
confidence-interval checks, but its late-half residual failed the fixed
absolute gate at both `t=60` and `t=80` (`0.01527` and `0.01013`,
respectively).  The latter is close to, but still above, `0.01`; the threshold
is not relaxed.

Before inspecting any further seeds, an independent remedial series is frozen
at `t=120` for both `N_c=512` (seed base `89600`) and `N_c=1024` (seed base
`89700`), with all other production parameters unchanged.  Acceptance at
`n=40` now requires the original gates at `t=120`, time convergence over the
recorded horizons through `t=120`, and population convergence at `t=120`.
The failed `t=60/80` series remains visible and is not pooled into the remedial
series.

## Post-summary n=30 population amendment and outcome

The first frozen cross-chain summary exposed failed individual-member
population convergence at `n=30`.  Before new data, the one allowed remedial
comparison was fixed at `N_c=1024 -> 2048`, with unchanged settings and gates;
see `REMEDIAL_AMENDMENT_2026-08-29.md`.  The final paired residual, support,
and time gates pass, but the `k=0.4` member change is about `2.36` combined SE
and therefore fails the unchanged population gate.  Per the stopping rule, no
additional population or tolerance is selected and `n=30` remains unresolved.

## Endpoint numerical-control amendments and outcomes

The initial `n=10`, `N_c=512` control baseline missed the frozen high-tilt ESS
floor, and selection interval 4 was unsupported.  Before replacement data,
`NUMERICAL_CONTROL_AMENDMENT_2026-08-29.md` fixed supported `N_c=1024`
timestep-halving and selection-interval-1 controls.  Both pass.

At `n=40`, timestep halving passes.  Selection interval 4 loses high-tilt
support.  Before replacement data,
`NUMERICAL_CONTROL_AMENDMENT_N40_2026-08-29.md` fixed one supported
selection-interval-1 comparison and prohibited further selection.  It fails
the unchanged individual-member/two-SE and late-half absolute gates.  The
core `n=40` paired result is therefore retained as positive but incomplete
evidence, not promoted to a fully controlled result.
