# Frozen protocol: one-bath open-chain stationary-profile diagnostic

Freeze date: 2026-09-11 (America/New_York).

## Boundary definition and implementation lineage

The open chain is defined by bath contributions, not a historical case
number.  The left endpoint uses the canonical BC1 bath,

`b_I,left = 2 gamma [2 T_left - (2 M I_1 - I_1^2 + 2 I_1 I_2 cos(delta_1))]`,

`b_phi,left = gamma [2 I_2 sin(delta_1)]`,

`sigma_I,left = 2 sqrt(2 gamma T_left I_1)`, and
`sigma_phi,left = sqrt(2 gamma T_left/I_1)`.

At the right endpoint the complete bath increment is absent:

`b_I,right=b_phi,right=sigma_I,right=sigma_phi,right=0`.

Hamiltonian drift at the terminal site remains.  The numerical implementation
is a prospective case 4 port onto the source of the completed 72-run profile
matrix, `boundary_profiles_2026-09-05/src/NLS_boundary_profiles.cpp`, whose
SHA-256 is
`40d13c7547d9c147ca8241ed701e1a1e85540fc0c21db48c4a63a943eccfe3ec`.
It retains that matrix's corrected noise, single-precision SIMD state,
adaptive Euler--Maruyama step `max(1e-5,min(5e-4,1/max_drift))`, polynomial
trigonometry, and positive-action projection.  Projection and invalid counts
are reported.  Because the completed controlled-integrator audit established
that this inherited scheme creates spurious bond-sine means, sine profiles are
recorded for completeness but are not quantitative evidence.

Prospective open-chain source SHA-256:
`dfee7ed0b416a30f2b27d919c524dfaf37e91ca44312ffbb9d68b7fc1cd36110`.
Compiled binary SHA-256:
`09d801e71b7f3bd13cea56aba3d432b471dbf16a726a9866894148a91397575d`.
The remote freeze commit is written to `FREEZE_COMMIT.txt` immediately after
the prospective source, protocol, scripts, and binary are pushed and before
production begins; the launcher copies it into the production manifest.

## Frozen run matrix

- `n = 25, 50, 100`; `gamma=0.1`; maximum `dt=5e-4`.
- Left-bath temperatures `T_left=10` and `T_left=6`; no right temperature.
- Burn-in `2000, 8000, 32000` at `n=25,50,100`.
- Measurement duration `2000`.
- Two independent seeds per condition, 16 SIMD batches x 16 lanes = 256
  trajectories per replicate and 512 trajectories after merging.
- Production is sequential with two OpenMP threads so it can occupy the two
  physical cores not used by the concurrently running FT audit.

Every output records the full boundary equations, source/binary identity,
seed, thread count, burn-in, measurement time, projection count, and sample
count.  Four burn-in snapshots and four post-burn state snapshots are fixed at
quarters of their intervals.  Four separate quarter-averaged profiles, rather
than only cumulative profiles, are saved.

## Predeclared stationarity diagnostics

Stationarity is evaluated before any scaling or collapse claim.  For the last
two separate measurement quarters (Q3 and Q4), define the symmetric relative
change as `2 |Q4-Q3|/(|Q4|+|Q3|)`.

A logical condition is called stationary only if all of the following hold:

1. the paired Q4-Q3 total-action difference has a 95% CI containing zero and
   relative change at most 5%;
2. the paired Q4-Q3 free-end-action difference has a 95% CI containing zero
   and relative change at most 5%;
3. across the interior `0.10 <= j/(n-1) <= 0.90`, the median sitewise relative
   Q3--Q4 change is at most 5% and the maximum is at most 10%.

All burn-in and measurement snapshots, total action, endpoint action, quarter
differences, z scores, and failures are reported.  A nonstationary profile is
retained as a time-dependent result; its scaling summaries are descriptive and
cannot be called stationary.

## Frozen profile analysis

Replicates are merged from their independent trajectory-level means and
variances.  Raw and `sqrt(n)`-rescaled profiles use pointwise 95% normal bands.
For each left temperature:

- fit terminal and midpoint action to `C n^alpha` using the three chain
  lengths; quote the one-residual-degree-of-freedom limitation;
- interpolate profiles linearly to the fixed 101-point grid
  `x=0,0.01,...,1` and report median/maximum relative spread over the full grid
  and the fixed interior `0.10<=x<=0.90`, for raw and rescaled profiles;
- define penetration depth prospectively as the first site where the merged
  mean action is at most 50% of the driven-end value.  Report the integer site,
  interpolated crossing, and its fraction of chain length.  If there is no
  crossing, report `not reached` without extrapolation.

No exponent, threshold, burn-in, fit window, or stationarity rule may be
changed after production output exists.
