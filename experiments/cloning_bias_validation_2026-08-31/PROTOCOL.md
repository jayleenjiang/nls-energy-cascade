# Prospective protocol: controlled-cloning population-bias validation

**Frozen before new simulation output:** 2026-08-31.

## Purpose

This experiment addresses three previously missing safeguards for the
controlled population-dynamics estimator of

```text
psi_n(k,t) = t^{-1} log E_NESS exp[-k Sigma_R(t)],
Sigma_R = -(1/T_R-1/T_L) Q_R = -0.4 Q_R,
```

at `(T_L,T_R,gamma)=(10,2,0.1)`.  It tests the estimator on the actual
interacting `n=2` dynamics, measures its finite-population movement, and
replaces single-population GC residuals by `1/N_c -> 0` intercepts.

No result from this experiment may be described as a proof.  If a reference,
support, plateau, or extrapolated-symmetry gate fails, the corresponding
long-chain GC statement is **not numerically established**.

## Observable matching for the n=2 known-answer test

The existing accepted short-time `n=2` result concerns total entropy

```text
Sigma_total = Sigma_medium + log rho_NESS(X_0) - log rho_NESS(X_t).
```

The controlled-cloning program used for the long chains instead weights the
right-bath gauge observable `Sigma_R=-0.4 Q_R`.  These two observables differ
by a random endpoint term at finite `t`, so their numerical SCGF values must
not be subtracted.  The accepted total-entropy values will be reproduced and
reported as a separate control, but the primary known-answer gaps use an
independent direct estimate of the **same** `Sigma_R` observable.

The new direct `n=2` reference is fixed prospectively at `t=0.1`, with 1024
independent streams and 1024 non-overlapping blocks per stream after burn-in.
The primary cloning validation uses the same horizon and observable at
`k=0.3,0.5,0.7`.  At this horizon the estimator has one fixed-population
selection event; it therefore validates the controlled transition kernel,
likelihood ratio, log-normalizer, and finite-`N_c` log-mean bias on the real
model.  It does not by itself validate long-time genealogy.

## Frozen simulator and common parameters

- source: `flux/NLS_entropy_cloning.cpp`;
- source SHA-256 before build:
  `bc3c27bf62a45aa879c4a8fd3e4d70fe7bdf1e3bcb0d06c0318826e61e615b6d`;
- Git source blob: `0f1160ec03f480582639c2d23f7043b1f6a39260`;
- mode: `controlled` (exact finite-step Gaussian likelihood correction);
- temperatures: `10,2`;
- burn-in: `500`;
- timestep: `5e-4`;
- gauge shift: `0.1`, hence heat coefficients `(0,-0.4)`;
- control scale: `0.5`;
- fixed resampling at every selection event;
- four independent seeds at every `(n,k,N_c)` cell;
- no seed is reused between cells.

The release executable is built once into this experiment directory, tested,
hashed, and used unchanged for every new run.  Its hash and compiler command
are recorded in `BUILD_PROVENANCE.md` before the first simulation starts.

## Frozen run matrix

### A. Physical-model known answer (`n=2`)

- direct reference: 1,048,576 `t=0.1` blocks;
- cloning tilts: `0.3,0.5,0.7`;
- populations: `512,1024,2048,4096`;
- horizon/selection interval: `0.1/0.1`.

### B. Long-chain population extrapolation

| n | tilt pair | horizon | selection | populations |
|---:|:---:|---:|---:|:---|
| 10 | `(0.3,0.7)` | 60 | 2 | 512,1024,2048,4096 |
| 20 | `(0.4,0.6)` | 60 | 2 | 512,1024,2048,4096 |
| 30 | `(0.4,0.6)` | 60 | 2 | 512,1024,2048,4096 |
| 40 | `(0.4,0.6)` | 120 | 2 | 512,1024,2048,4096 |

The `n=10,k=0.3` timeseries checkpoint at `t=20` is also compared with the
already frozen direct value

```text
-0.182236735448602,
95% stream-bootstrap CI [-0.18269705679294015,-0.18170401300354538].
```

This checkpoint comes from the same `t=60` runs; no value is selected from the
timeseries after inspection.

All primary seeds and exact command arguments are generated in
`RUN_MATRIX.csv`.  Conditional `N_c=8192` seeds are also frozen there before
execution.

## Support and integrity gates

Every accepted run must have a summary and timeseries, finite numeric fields,
the exact frozen source/binary identifiers, zero midpoint failures, the
expected final time, and minimum selection-weight ESS at least `0.1 N_c`.
Missing or failed runs are never silently retried or replaced with another
seed.

## Primary `1/N_c` model and uncertainty

For each fixed `(n,k,t)`, all per-seed estimates are fit without averaging
away the raw runs:

```text
psi_r(N_c) = a + b/N_c + epsilon_r.
```

The `N_c=infinity` estimate is `a`.  The primary covariance is the HC3 robust
OLS covariance across all independent seed estimates.  A 95% Student interval
uses `df = number_of_runs - 2`.  Complementary tilts use disjoint seeds, so
the variance of `a(k)-a(1-k)` is the sum of the two intercept variances; the
reported interval uses a Welch--Satterthwaite degree of freedom.  A fixed-seed
hierarchical bootstrap is a sensitivity analysis only.

No built-in correction is claimed: the simulator reports the usual
finite-population log-normalizer.  The external intercept is the prospective
finite-population correction.

## Population plateau rule and conditional `N_c=8192`

The first decision uses only each individual `psi(k)` series, never the GC
residual and never closeness to a reference.  The frozen `n=10,k=0.3,t=20`
checkpoint is included as an additional individual-member series.  A member
passes at 4096 only if:

1. all four population levels and four seeds per level pass support;
2. `|mean(4096)-mean(2048)| <= max(0.002, 2 combined SE)`;
3. the fitted correction at 4096 obeys `|b|/4096 <= 0.002`;
4. the HC3 intercept SE is at most `0.0025`.

If any member at a given `n` fails, **all members fixed for that n** are run at
`N_c=8192` using the already listed seeds.  Thus an additional population is
never selected because it makes a pair look more symmetric.  After 8192 there
is no further adaptive population choice: failure to plateau is reported as
unresolved.

## Final decision rules

The estimator validation requires all of the following:

1. at each `n=2` tilt, the extrapolated same-observable gap to the direct
   `Sigma_R` reference has a 95% CI containing zero and absolute point gap at
   most `0.005`;
2. at `n=10,t=20,k=0.3`, the extrapolated gap to the frozen direct reference
   has a 95% CI containing zero and absolute point gap at most `0.003`;
3. every long-chain member passes the final support and plateau audit;
4. each extrapolated complementary-pair residual has a 95% CI containing zero,
   absolute point residual at most `0.005`, and CI half-width at most `0.005`.

If rule 1 or 2 fails, the cloning estimator is not validated at the precision
needed for the reported residuals, and no long-chain GC symmetry claim stands.
If rules 3 or 4 fail for a chain length, GC symmetry is not numerically
established for that length.  A negative or unresolved verdict is final under
this protocol.

## Execution order

1. archive/build/hash/self-test the executable;
2. generate the new direct `n=2` reference;
3. run all primary `n=2` and long-chain population cells;
4. run the member-only plateau audit;
5. run the predeclared 8192 cells only where requested by that audit;
6. run final analysis, integrity checks, and report generation.

The AC-wait runner checks for an existing matching process and never launches
a duplicate.  Loss of AC power never kills an active pair; it delays the next
pair.
