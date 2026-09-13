# Prospective protocol: long-chain Gallavotti--Cohen Phase II

**Frozen:** 2026-08-29, before generating any Phase-II simulation output.

## Research question and scope

For the projection-free Cartesian NLS chain at
`(T_L,T_R,gamma)=(10,2,0.1)`, determine whether the resolved entropy-current
SCGF pair

```text
psi_n(k) = psi_n(1-k),  k=0.4,
Sigma_R = -(1/T_R-1/T_L) Q_R = -0.4 Q_R
```

passes the missing population, timestep, and selection-interval convergence
checks at `n=20,30,40`.  Passing supports numerical consistency over this
tested pair.  It is not a proof, an all-tilt statement, or an infinite-chain
limit.

This is a new prospective validation phase.  It does not retroactively change
the frozen Phase-I outcomes or thresholds.

## Existing evidence reused without refitting

- Direct two-tail data: `experiments/entropy_ft_2026-08-26/production/`.
  Each `n=10,20,30,40` file contains 1,000,064 non-overlapping `t=20` blocks,
  aggregated to `t=20,40,...,200`.  The accepted audit is
  `production/final_v4/` and the compact results are in `curated_results/`.
- Phase-I cloning baselines:
  - `n=20`: `N_c=1024`, `t=60`, `dt=5e-4`, selection interval 2;
  - `n=30`: `N_c=2048`, `t=60`, `dt=5e-4`, selection interval 2;
  - `n=40`: `N_c=1024`, `t=120`, `dt=5e-4`, selection interval 2.
- The source is `flux/NLS_entropy_cloning.cpp` at the Git revision recorded
  before execution.  The exact Gaussian controlled-kernel and cloning
  self-tests must pass before production.

## Fixed Phase-II matrix

Every cell uses four new independent seeds for each of `k=0.4` and `k=0.6`.
The two tilt members may run concurrently; seed groups are never reused.

| cell | n | N_c | horizon | dt | selection | seeds | purpose |
|---|---:|---:|---:|---:|---:|---:|---|
| N20-DT | 20 | 1024 | 60 | 2.5e-4 | 2 | 91001--91008 | timestep halving |
| N20-SEL | 20 | 1024 | 60 | 5e-4 | 1 | 91101--91108 | selection interval |
| N30-POP | 30 | 4096 | 60 | 5e-4 | 2 | 92001--92008 | population doubling |
| N40-SEL | 40 | 1024 | 120 | 5e-4 | 1 | 93001--93008 | selection interval |
| N40-DT | 40 | 1024 | 120 | 2.5e-4 | 2 | 93101--93108 | timestep halving |
| N30-DT | 30 | 4096 | 60 | 2.5e-4 | 2 | 92101--92108 | conditional timestep control |
| N30-SEL | 30 | 4096 | 60 | 5e-4 | 1 | 92201--92208 | conditional selection control |

Common settings are burn-in 500, `gauge_shift=0.1`, `control_scale=0.5`, fixed
population resampling at every selection event, and five OpenMP threads per
process.  No control scale, tilt pair, horizon, population, or tolerance is
selected from Phase-II results.

## Hard numerical and support gates

Every accepted group must satisfy all of the following:

1. at least four independent seeds per tilt member;
2. zero nonfinite trajectories and zero midpoint failures;
3. minimum selection-weight ESS at least `0.1 N_c` in every run;
4. finite full-window and late-half SCGF estimates;
5. the recorded heat/entropy gauge identity and first-law checks retain the
   existing numerical tolerance.

## Symmetry and convergence gates

The Phase-I thresholds are retained unchanged:

- the 95% independent-run interval for
  `D_n(k)=psi_n(k)-psi_n(1-k)` contains zero;
- `|D_n(k)| <= 0.01`, with the same conclusion for the late-half slope;
- changing population, timestep, or selection interval changes each member,
  the paired residual, and the late-half paired residual by no more than
  `0.01` and two combined independent-run standard errors;
- the final pair retains the Phase-I observation-time convergence gate.

The direct-sampling probability ratios are descriptive finite-time evidence.
Plus-four smoothing and Gaussian tail extrapolation are shown only as
diagnostics and are never used to manufacture a resolved FT tail.

## Conditional stopping rule

Stage I runs N20-DT, N20-SEL, N30-POP, N40-SEL, and N40-DT.  N30-DT and
N30-SEL are run only if N30-POP passes both individual-member population gates,
the paired population gate, support, time convergence, and the final paired
symmetry gate.  If N30-POP fails, `n=30` remains unresolved and its controls
are not run.

No failed Phase-II cell is followed by another data-dependent population,
selection interval, control scale, tilt pair, or relaxed tolerance.  Such a
change would require a separately dated protocol and would remain a new
experiment rather than a repair of this one.

## Stronger secondary audit

For interpretation, report the normalized residual

```text
R_n = |psi(k)-psi(1-k)| / ((|psi(k)|+|psi(1-k)|)/2)
```

and its confidence interval.  This is a secondary effect-size diagnostic, not
a replacement for the frozen primary gates.  A non-significant residual alone
is not described as proof of equality.
