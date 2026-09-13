# Mentor update: finite-time current fluctuations

## One-sentence result

The requested million-sample, two-tail data collection is complete and fully
audited.  It is production-quality for the reported finite-time observables
but remains provisional as an FT test: the action-current tails show a
well-resolved finite-time
large-deviation scaling `P ~ exp[-t I(A)]`, but neither the medium-entropy nor
the bath-heat symmetry has reached its simple asymptotic fluctuation-theorem
reference in the directly sampled windows, so the correct conclusion is
“not verified at accessible times,” not “the FT is proved” or “the FT is
violated.”

## What was simulated

- Projection-free Cartesian dynamics `c_j=x_j+i y_j`; no action floor.
- `T_left=10`, `T_right=2`, `gamma=0.1`, `dt=5e-4`, burn-in `B=500`.
- `n=10,20,30,40`; 128 independent streams per chain length.
- 1,000,064 non-overlapping base blocks of length `t=20` per chain length.
- Blocks were aggregated to `t=20,40,...,200` without inventing zero-count
  tails.
- Every block stores `Q_left`, `Q_right`, `Delta E`, medium entropy,
  action current, stream ID, block ID, and the first-law residual.

The heat sign is positive into the system and

`Sigma_medium = -Q_left/T_left - Q_right/T_right`.

## Numerical validity

- Four raw files contain exactly 1,000,064 finite, correctly ordered rows.
- Entropy and first-law columns recompute exactly from the raw heat columns.
- Midpoint failures: zero.
- Balance RMS rate decreases from `9.71e-6` at `n=10` to `3.05e-6` at
  `n=40`; halving `dt` reduces it by approximately a factor of four.
- Equal-temperature control has zero current within uncertainty; bath swap
  reverses the currents while mean entropy remains positive.
- Maximum first/last-quarter stationarity diagnostic is only `|z|=1.34` over
  all observables and requested windows.
- The final independent analysis audit performed 75,617 checks with zero
  errors.

## Mean currents

| n | action current | heat current | medium entropy rate |
|---:|---:|---:|---:|
| 10 | 0.39809 | 3.21212 | 1.28485 |
| 20 | 0.11979 | 1.40691 | 0.56275 |
| 30 | 0.05426 | 0.79595 | 0.31839 |
| 40 | 0.03052 | 0.52294 | 0.20918 |

The joint Cartesian sampler gives

`<J_M(n)> = 28.889 n^(-1.8485)` (`R^2=0.9987`),

which independently agrees with the established canonical result
`28.75 n^(-1.850)`.

## Dependence on threshold A

At fixed averaging time, `log P[J_t >= A]` and `log P[J_t <= -A]` are curved
functions of `A`; a single exponential-in-`A` law does not describe the whole
tail.  The central 1--99% region becomes increasingly Gaussian-looking with
time, but the `t=20` far tails are heavier than the Gaussian fixed by the full
sample variance.  A joint two-tail Gaussian fit requires a width about
13%, 20%, 26%, and 28% larger than the full-sample width for
`n=10,20,30,40`, respectively.  The Gaussian fit is therefore a descriptive
benchmark, not an FT test.

## Dependence on averaging time t

For raw-count-qualified rare-tail thresholds we fitted

`log P_t = intercept - t I_+(A)` and
`log P_t = intercept - t I_-(A)`.

All resolved fits have `R^2 >= 0.98`:

| n | positive-tail A range | fits | negative-tail A range | fits |
|---:|---:|---:|---:|---:|
| 10 | [0.46, 0.81] | 36 | unresolved | 0 |
| 20 | [0.15, 0.33] | 19 | [0.01, 0.06] | 6 |
| 30 | [0.08, 0.20] | 13 | [0.01, 0.08] | 8 |
| 40 | [0.05, 0.15] | 11 | [0.01, 0.08] | 8 |

The fitted rate proxies increase with threshold magnitude.  This is strong
finite-time large-deviation evidence.  The `n=10` negative tail becomes too
rare after short times for a four-time-point direct-sampling fit; plus-four
smoothing is not used to claim otherwise.

## Fluctuation-symmetry result

For medium entropy, the simple long-time reference is

`t^(-1) log[p_t(a)/p_t(-a)] = a`.

The primary fitted slopes are below one.  They generally rise with averaging
time where the negative tail remains resolved—for example, at `n=30` they are
0.247, 0.316, and 0.393 for `t=20,40,60`—but do not reach one before direct
sampling loses two-sided resolution.  The heat-current slopes are likewise
below their reference `Delta beta=0.4`.

This does not establish a violation.  An exact finite-time detailed FT applies
to total trajectory entropy

`Delta s_total = Sigma_medium + Delta s_system`.

The NESS system-entropy endpoint term is not reconstructed here, and boundary
energy terms can strongly affect heat tails.  The allowed statement is:

> The standard thermodynamic fluctuation symmetry was not verified in the
> directly sampled finite-time window; the observed slopes move in the
> expected direction for several chain lengths but remain boundary-term and
> rare-event limited.

The action-current symmetry slope has no thermodynamic target.  It crosses
values near one for some `n,t`, but it is visibly fit-window dependent, so this
cannot be presented as a fluctuation theorem.

## Heat--action coupling

Correlation strengthens with averaging time, but tight coupling is not
established.  At `t=200`, Pearson correlations are 0.934, 0.896, 0.803, and
0.683 for `n=10,20,30,40`; the unexplained heat-variance fractions remain
0.128, 0.197, 0.355, and 0.533.  Action-current fluctuations therefore cannot
be substituted for entropy production, especially for longer chains.

## What to tell the mentor

“上次您让我检查 flux 两侧尾部以及它对时间窗口的依赖。现在
`n=10,20,30,40` 每个都有 1,000,064 个 `t=20` 样本，并聚合到了
`t=200`。固定阈值时，resolved tails 的 `log P` 对 `t` 基本线性，说明
存在清楚的 finite-time large-deviation scaling；rate function 随 `A`
增大。中心看起来接近 Gaussian，但 far tails 比 full-sample Gaussian
更重。对于 fluctuation theorem，medium entropy 和 heat 的 symmetry
slope 随 `t` 有向理论值移动的趋势，但在负尾消失前仍没有达到 1 和
0.4，所以目前只能说 FT 在可采样窗口内没有被验证。Action current
不能直接当 entropy production，因为 joint data 显示二者并非 tight
coupling。所有 raw、timestep/control 和独立复算 audit 都通过了。”
