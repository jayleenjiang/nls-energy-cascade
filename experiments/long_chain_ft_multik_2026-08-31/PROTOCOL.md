# Prospective protocol: multi-tilt long-chain GC symmetry

**Frozen:** 2026-08-31, before any Phase-III simulation output was generated.

## Question and claim boundary

For the projection-free Cartesian boundary-driven NLS chain at

```text
(T_L,T_R,gamma) = (10,2,0.1),  burn-in = 500,
Sigma_R = -(1/T_R-1/T_L) Q_R = -0.4 Q_R,
```

test the Gallavotti--Cohen SCGF relation

```text
psi_n(k) = psi_n(1-k)
```

at several complementary tilt pairs for each of `n=10,20,30,40`.  Phase II
resolved one pair per length.  Phase III asks whether the same conclusion is
stable over a predeclared, support-resolved multi-`k` window.

A passing result is a controlled finite-size numerical consistency statement.
It is not an exact equality, an all-tilt result, a proof uniform in `n`, or an
`n -> infinity` theorem.

## Frozen implementation

- source: `flux/NLS_entropy_cloning.cpp`;
- executable: `flux/entropy_cloning_v2`;
- mode: exact controlled Gaussian-kernel cloning (`controlled`);
- common settings: `dt=5e-4`, selection interval `2`, gauge shift `0.1`,
  control scale `0.5`, five OpenMP threads per process;
- each complementary pair is run concurrently, but successive seed pairs are
  run serially;
- every completed run must have both a nonempty summary and timeseries.

The source and script hashes are recorded before the pilot starts.  Existing
Phase-II files are read-only baselines and are never overwritten.

## Stage A: support-only pilot

The candidate pairs are

```text
(0.25,0.75), (0.35,0.65), (0.45,0.55).
```

Every candidate is tested at each chain length with `N_c=512`, horizon `20`,
and one independent seed per tilt member.  The exact seeds are in
`PILOT_MATRIX.csv`.

The pilot grid is selected **without reading the SCGF estimate or pair
residual**.  A tilt member is support-eligible only if:

1. the summary and timeseries are present and finite;
2. `midpoint_failures = 0`;
3. minimum weight ESS is at least `0.1 N_c`;
4. minimum unique roots is at least `32`;
5. minimum root-count ESS and root-weight ESS are each at least `16`.

A pair is eligible only if both members pass.  For each `n`, the deterministic
selection rule chooses the outermost and innermost eligible pairs.  If fewer
than two pairs are eligible, that chain length is declared blocked and no
Phase-III production is started for any length.  No symmetry value may be
used to rescue, replace, or reorder a candidate.

The resulting `FROZEN_GRID.csv` is the second-stage prospective production
matrix.  It is created once and must not be edited after production begins.

## Stage B: production

Four new independent seeds are run for each member of each selected pair.
The final settings reuse the converged Phase-II populations and horizons:

| `n` | `N_c` | horizon | existing Phase-II pair |
|---:|---:|---:|:---:|
| 10 | 2048 | 60  | `(0.3,0.7)` |
| 20 | 1024 | 60  | `(0.4,0.6)` |
| 30 | 4096 | 60  | `(0.4,0.6)` |
| 40 | 1024 | 120 | `(0.4,0.6)` |

The existing pair is reused without refitting.  Together with the two new
pairs, each chain length therefore has three resolved complementary pairs.

## Stage C: controls

For the outermost newly selected pair at every length, four seeds per member
are run at half the production population.  At `n=40`, the same outer pair is
also run with:

- timestep halving, `dt=2.5e-4`, at the production population;
- selection interval `1`, at the production population and baseline timestep.

These controls are run irrespective of the observed symmetry residual, as
long as the production summaries exist and pass basic numerical support.

## Frozen gates

Every accepted production group must have four independent seeds, zero
midpoint failures, finite full-window and late-half estimates, and minimum
weight ESS at least `0.1 N_c` in every run.  At every selected pair:

1. the independent-run 95% interval for
   `D_n(k)=psi_n(k)-psi_n(1-k)` contains zero;
2. `|D_n(k)| <= 0.01` for both full-window and late-half estimates;
3. the final-vs-penultimate observation-time change is at most `0.01` and two
   combined independent-run standard errors;
4. at the outer pair, both tilt members and the paired residual pass the same
   population-convergence tolerance;
5. the `n=40` outer pair passes both timestep and selection-interval controls
   for each member and for the paired residual.

The already frozen Phase-II gates remain the acceptance criteria for the
reused baseline pair.  A failed Phase-III cell is reported as failed; no
data-dependent change of `k`, population, horizon, control scale, numerical
tolerance, or seed count is permitted.

## Outputs

The final analysis will contain the complete run matrix, per-pair residuals
with confidence intervals, support and convergence gates, a multi-`k` SCGF
plot, a residual plot, an integrity report, and a claim-bounded final verdict.

