# Frozen protocol: direct finite-time entropy test at n=3

Frozen before production output is inspected.  The direct-sampling stage is
the decision gate; no cloning run belongs to this protocol.

## Scientific question

At `n=3`, `(T_left,T_right)=(10,2)`, and `gamma=0.1`, determine whether
ordinary direct sampling resolves the negative medium-entropy tail at block
horizons `t=20,40,...,200`.  If it does, report the directly observable
medium-entropy symmetry diagnostic.  A finite-time total-entropy test is made
only if a stationary density for these exact parameters passes an independent
density audit.

## Frozen dynamics and sampling parameters

- Dynamics: projection-free Cartesian Gibbs-preserving NLS dynamics.
- Integrator: alternating left/right bath order; implicit midpoint
  Hamiltonian step and separate stochastic bath steps, identical to the
  `n=2` endpoint and `n=10,...,40` heat samplers.
- `n=3`, `T_left=10`, `T_right=2`, `gamma=0.1`.
- `dt=0.0005`, burn-in `500`, base block time `20`.
- `8` SIMD batches x `16` lanes = `128` independent streams.
- `7813` non-overlapping base blocks per stream, hence `1,000,064` base
  blocks.
- Production seed: `2026083133`.
- Recorded columns: stream and block identifiers, `Q_left`, `Q_right`,
  `Delta E`, `Sigma_m=-Q_left/10-Q_right/2`, entropy rate, middle-bond action
  current, first-law residual, and the reduced five-dimensional state
  `(I1,I2,I3,theta1,theta3)` at both endpoints.
- Reduced angles are `theta1=2(phi1-phi2)` and
  `theta3=2(phi3-phi2)`, wrapped to `[-pi,pi)`.

The first `floor(7813/m)*m` base blocks in each stream are grouped without
overlap for horizon `t=20m`; leftover blocks are discarded only for that
horizon.  Streams are never concatenated.

## Frozen feasibility rule

At each horizon, direct two-sided support is called *resolved* only when both
conditions hold:

1. at least `1000` raw blocks have `Sigma_m < 0`; and
2. at least `8` symmetric histogram-bin pairs contain at least `20` raw
   samples on each side.

Symmetric bins are centered at zero.  Their width is the unrounded
Freedman--Diaconis width computed from all blocks at that horizon.  All
qualifying pairs are reported; no interval is selected for closeness to the
reference slope.  The `t=20` rule decides whether direct sampling can support
a finite-time two-tail analysis.  Longer horizons are reported individually
and may lose support.

Negative-tail probability intervals use `2000` independent-stream bootstrap
replicates with analysis seed `2026083193`.

## Frozen medium-entropy analysis

At every horizon with at least two qualifying symmetric bins, the descriptive
medium-entropy relation is fitted to

```
log[N(+a)/N(-a)] = intercept + slope * a.
```

The representative `a` is the center of the symmetric FD bin pair.  The
primary fit is weighted least squares with free intercept and inverse
log-count-ratio variance weight
`1/(1/N_plus+1/N_minus)`.  Every full-sample bin pair with at least `20`
counts on both sides enters; this set is fixed for the stream bootstrap and no
contiguous subwindow is selected.  Report slope, intercept, weighted `R^2`,
all raw paired counts, and `2000` stream-bootstrap percentile intervals.

For the medium-only integral diagnostic, report
`log mean exp(-Sigma_m)`, its stream-bootstrap interval, exponential-weight
effective sample size, and the largest single-sample share of the exponential
sum.  These tail diagnostics determine whether the numerical IFT average is
actually resolved.  The reference values `slope=1` and `log IFT=0` belong to
total entropy; medium-only discrepancies are reported, not called FT
violations.

## First-law and integrity gates

- Exactly `1,000,064` finite, ordered base-block rows.
- Exact recomputation of `Sigma_m=-Q_left/10-Q_right/2` to floating-point
  output tolerance.
- Exact recomputation of the saved first-law residual
  `Q_left+Q_right-Delta E` to floating-point output tolerance.
- Zero midpoint failures.
- Consecutive block endpoints agree within output roundoff for every stream.
- All actions are nonnegative and all saved reduced angles lie in
  `[-pi,pi)`.
- The residual distribution is reported without deleting outliers.

## Total-entropy branch

The saved Section-4 nonequilibrium NN is not admissible for this production:
it was trained at `(T1,T3)=(2,8)`, not `(10,2)`, and its documented density
and Fokker--Planck errors are too large to treat it as an exact endpoint
density.  This incompatibility was found before production and is recorded in
`MODEL_COMPATIBILITY_AUDIT.md`.

Therefore:

- the existing Section-4 NN will **not** be evaluated on the `10,2` blocks;
- no total-entropy slope or IFT number will be manufactured from it;
- a total-entropy branch can start only after an independently trained
  `10,2` density model is frozen and validated on held-out equilibrium and
  nonequilibrium data, including endpoint log-density-ratio error;
- if no such model passes, the final report contains only the medium-entropy
  result and negative-tail resolution, as required by the honesty clause.

## Claim boundary

Even a successful direct two-tail result is a numerical finite-time result at
the stated parameters.  Medium entropy alone is not the exact finite-time
total entropy.  Failure of direct support is not evidence against the FT.  No
cloning result may be substituted silently.
