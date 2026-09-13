# Frozen protocol: n=10 heat-flux affinity sweep

Protocol freeze date: 2026-09-03 (America/New_York)

## Question

At fixed arithmetic mean bath temperature 6, test whether the long-time heat
fluctuation-relation slope approaches

`Delta beta = 1/T_R - 1/T_L`

in cases where the negative heat tail remains directly sampled.  In parallel,
measure the approach of the integrated heat distribution toward a Gaussian and
distinguish the Gaussian bulk slope `2*mean(Q)/Var(Q)` from the directly sampled
two-tail slope.

This is a numerical finite-time/finite-sample test.  A failed or unresolved gate
is retained as such; it is not replaced by smoothed, fitted, plus-four, or
extrapolated tail mass.

## Frozen dynamics and source

- Chain length: `n=10`.
- Cartesian, projection-free sampler: `NLS_entropy_ft.cpp`.
- Exact source commit used by the accepted `(T_L,T_R)=(10,2)` production:
  `1905cf4e606a4a7f4dd8930caa64bd4cc861e9d4`.
- Frozen source copy: `source/NLS_entropy_ft_1905cf.cpp`.
- Frozen source SHA-256:
  `98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3`.
- Model version emitted by the source: `gibbs-cartesian-entropy-ft-v1`.
- `gamma=0.1` (compile-time constant in this exact source).
- Integrator step: `dt=5e-4`.
- Burn-in: `500`.
- Base block duration: `20`.
- Bond: `5`.
- No source or integrator change is permitted after production begins.

The current repository source has a different hash because later commits added
n=3 endpoint output.  That newer source is not used for this sweep.

## Frozen cases, sample size, and seeds

Each new case uses 8 SIMD batches x 16 lanes = 128 independent streams and
7,813 non-overlapping base blocks per stream, for exactly 1,000,064 base blocks.
Two OpenMP threads are allocated per case.  The four cases are launched
concurrently only on AC power.

| case | T_L | T_R | Delta beta | seed |
|---|---:|---:|---:|---:|
| `dbeta_0p027972` | 6.5 | 5.5 | 0.02797202797202797 | 2026090401 |
| `dbeta_0p057143` | 7 | 5 | 0.05714285714285716 | 2026090402 |
| `dbeta_0p125000` | 8 | 4 | 0.125 | 2026090403 |
| `dbeta_0p222222` | 9 | 3 | 0.2222222222222222 | 2026090404 |

The exact existing `(10,2), n=10` production is reused from
`/Users/jayleenjiang/Downloads/n10_blocks.csv`; its SHA-256 is
`a23806e82f5514a9c3375d10a6644946b6b57d0efe8452b3bbf397b6230f9929`.

## Equilibrium availability audit

No `(6,6), n=10, 1,000,064-block` data set was found locally or in the tracked
experiment manifests.  The available equal-temperature files are:

- `(6,6), n=20`, 1,280 blocks: a timestep/control sample, not the requested
  chain length or production size;
- `(6,6), n=3`, 1,000,064 blocks: a known-answer heat audit, not `n=10`.

Neither is scientifically interchangeable with the requested `n=10` CLT row.
Because the instruction explicitly says not to rerun equilibrium, this protocol
does not launch a fifth case.  The `Delta beta=0` crossover row is reported as
`UNAVAILABLE` rather than silently substituting a different chain length.  A
future `n=10` equilibrium production requires explicit authorization and a new
protocol amendment made before that run.

## Frozen observable and aggregation

For each base block,

`Q = (Q_left - Q_right)/2`.

Non-overlapping windows `t=20,40,80,160,320,640` are formed only by summing
consecutive base blocks within the same stream.  If 7,813 is not divisible by
the number of base blocks per window, the trailing remainder is discarded
separately in each stream.  Streams are never joined across a window.

## Frozen per-window statistics

For every available `(case,t)` report the raw window count, mean, sample
standard deviation (`ddof=1`), moment skewness, moment excess kurtosis, and
`n_neg = #{Q<0}`.  Also report

- `a_Gauss = 2*mean(Q)/Var(Q)`;
- `gaussFT = Var(Q)/(2*mean(Q)/Delta beta)` for nonzero `Delta beta`.

The Gaussian quantities are bulk diagnostics, not substitutes for directly
sampled two-tail evidence.

## Frozen matched-bin two-tail fit

The rule is identical for every case and time:

1. Set `dx = sample_std(Q)/20`.
2. Use bins centred at integer multiples `k*dx`.
3. Mark a bin reliable only when its raw count is at least 10.
4. Among contiguous reliable-bin blocks that strictly straddle zero, select
   the block with the largest number of bins.  If tied, select the block with
   the smaller left index.  No block that fails to straddle zero is eligible.
5. Reflect the selected block to its largest symmetric range around zero.
6. Retain a positive/negative bin pair only if both raw counts are at least 10.
7. Require at least three retained positive/negative pairs; otherwise mark the
   result `UNRESOLVED`.
8. Fit `log(c_plus/c_minus) = a_fit*Q + intercept` by weighted least squares,
   with variance estimate `1/c_plus + 1/c_minus` for each ratio point.

There is no plus-four correction, KDE, normal-tail extrapolation, or adaptive
per-case tuning in this fit.

The analytic WLS slope SE is reported.  A separate stream bootstrap uses 1,000
multinomial resamples of the 128 complete streams.  The full-sample `dx`, bin
edges, and accepted symmetric window remain fixed during bootstrap so the
bootstrap measures stream-level sampling uncertainty conditional on the frozen
fit construction; the number of resolved bootstrap replicates is reported.
Percentile 95% intervals require at least 800 resolved replicates.

## Frozen long-time extrapolation and plateau rule

If at least three full-sample `a_fit` values are resolved, fit

`a_fit(t) = a_inf + c/t`

by unweighted ordinary least squares using every resolved requested time.  The
95% interval for `a_inf` is obtained by a joint stream bootstrap: each bootstrap
replicate uses the same resampled streams at every time and the same
full-sample bins/windows.  The reported FT comparison is `a_inf/Delta beta`
with the propagated percentile interval.  No time point is dropped because it
makes the intercept less favourable.

Visible plateau is `YES` only if the two largest resolved times are consecutive
doublings and

`abs(a_last-a_previous)/abs(a_previous) < 0.05`.

Otherwise it is `NO`; if fewer than two appropriate times are resolved it is
`UNAVAILABLE`.

## Frozen crossover summary and figures

The summary contains the six requested affinities `0, 0.027972, 0.057143,
0.125, 0.222222, 0.4`.  The equilibrium row remains `UNAVAILABLE` under the
availability audit above.  For each nonzero case report:

- `a_inf/Delta beta`, or `UNRESOLVED`;
- `a_Gauss/Delta beta` at the largest time with a resolved two-tail fit;
- `gaussFT` at `t=640`;
- skewness and excess kurtosis at `t=160`;
- `n_neg` at `t=160`.

Figures are fixed as:

1. `a_inf` and the selected `a_Gauss` versus `Delta beta`, with `y=Delta beta`;
2. `gaussFT(t=640)` versus `Delta beta`.

## Frozen decision language

- If the 95% CI for `a_inf/Delta beta` excludes 1, the FT test for that case is
  recorded as `FAIL`.
- If it includes 1, the finite-sample extrapolation is recorded as
  `CONSISTENT_WITH_FT`, without claiming proof.
- If fewer than three times or fewer than 800 joint-bootstrap intercepts are
  resolved, the result is `UNRESOLVED`.
- A negative or null result is retained without post-hoc physical explanation.
