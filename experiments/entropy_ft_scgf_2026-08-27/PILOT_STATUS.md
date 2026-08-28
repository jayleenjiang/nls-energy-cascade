# SCGF / fluctuation-theorem pilot status

## Outcome

The numerical machinery is validated in the low-tilt region, but a
Gallavotti--Cohen fluctuation-theorem comparison is **not yet accepted**.  No
nontrivial pair `k` and `1-k` simultaneously has reliable rare-event support.

The target is the finite-time medium-entropy SCGF

`psi_t(k) = log E[exp(-k Sigma_medium)] / t`.

The asymptotic FT diagnostic would be `psi(k)=psi(1-k)`.  Equality at `k=1/2`
is algebraic and is not counted as evidence.

## Direct million-sample result

The direct estimator was applied to the audited `1,000,064` base blocks for
each of `n=10,20,30,40`.  All `14,762` independent audit checks passed.
However, exponential reweighting becomes dominated by a few samples.  At
`k=1`, the sample-weight ESS is only about `1--8` depending on chain length and
window; one sample can carry `28%--96%` of the total weight.  There is no
nontrivial `k <-> 1-k` pair for which both direct estimates pass the ESS and
stream-concentration gates.

## Cloning validation

- Gaussian exact-model self-test: pass.
- `k=0` recovery with four independent runs:
  - medium entropy rate `1.28906 +/- 0.00782`, versus direct `1.28485`;
  - action current `0.39897 +/- 0.00162`, versus direct `0.39809`.
- One-selection-interval identity between the naive path-weight and guided
  tilted-generator decompositions:
  - absolute SCGF difference `0.00290` at `k=0.1`;
  - absolute SCGF difference `0.01393` at `k=0.5`.
- At `n=10`, `t=20`, `k=0.3`, four-run guided estimates converge with clone
  population:

| clone population | mean SCGF | run SE | minimum unique roots | minimum root ESS |
|---:|---:|---:|---:|---:|
| 256 | -0.15685 | 0.00884 | 9 | 3.70 |
| 512 | -0.16985 | 0.01100 | 18 | 9.09 |
| 1024 | -0.16035 | 0.00513 | 40 | 20.87 |

The direct million-sample estimate is `-0.15784`.  The `N_c=1024` result
differs by `0.49` independent-run SE and passes the predeclared support gate.
The `N_c=512 -> 1024` change is `0.00950`, below both two combined SE
(`0.02429`) and the absolute `0.02` tolerance.

At `N_c=1024`, changing the selection interval from `Delta=0.5` to `0.25`
or `1.0` changes the mean SCGF by `0.00869` and `0.00452`, respectively.
Both differences pass the predeclared two-SE and absolute `0.01` criteria;
all three settings retain at least 35 initial roots.

Halving the timestep gives `-0.15866 +/- 0.00324`, compared with
`-0.16035 +/- 0.00513` at `dt=5e-4`.  The difference is `0.00169`, below two
combined SE (`0.01214`) and the absolute `0.01` criterion.  The fine run has
zero midpoint failures and retains at least 39 roots.

## Why the FT gate is still blocked

At `k>=0.5`, the tested populations lose genealogical diversity.  For the
standard `N_c=1024`, `t=20`, `k=0.7` attempt, three runs finish but each ends
with only one initial root; a fourth run is rejected by the midpoint numerical
gate before a summary is written.  The minimum within-selection weight ESS is
only `8.65`.  The corresponding direct estimate also fails its support gate.
These estimates are rejected rather than interpreted physically.

Extending the accepted low-tilt point to `t=40` gives
`psi_40(0.3)=-0.17607 +/- 0.00719` over four runs.  The endpoint populations
retain only `14--19` distinct initial roots and the minimum root ESS is `8.90`,
below the predeclared requirements of 32 roots and root ESS 16.  Longer
observation therefore worsens the ancestry problem even at `k=0.3`.

This failure is scientifically informative.  Medium entropy omits the NESS
system-entropy endpoint term.  In an unbounded state space, endpoint-energy
fluctuations can control exponential moments near `k=1`, so simply collecting
more ordinary samples does not guarantee an FT test.

## Small-chain total-entropy control

The separate `n=2` endpoint pilot generated `32,768` non-overlapping `t=20`
blocks from 64 streams for both equal-temperature `(6,6)` and driven `(10,2)`
conditions.  Both raw files pass independent formula, first-law, finiteness,
and midpoint audits; their first-law RMS rates are `1.24e-5` and `1.51e-5`.

The smoothing scale was selected on the first half of the equal-temperature
streams only.  On the held-out half, the learned stationary-density ratio has
support `0.9868`, slope `0.9296`, correlation `0.9675`, and RMSE `0.3603`
against the exact Gibbs ratio.  These ordinary-error gates pass, and using the
exact Gibbs endpoint gives `log mean exp(-Delta s_total)=-3.29e-7`.

The stricter learned-density IFT gate does not pass: its equilibrium value is
`0.0843` with stream-bootstrap interval `[0.0700,0.1025]`, which excludes zero.
For the driven data, the exponential-weight ESS is only `1.63`, and no
positive/negative histogram pair has 20 raw counts on both sides.  The learned
endpoint and driven symmetry are therefore blocked.  This is an estimator and
rare-event limitation, not evidence that the FT is violated.

## Allowed statement

> Direct and population-dynamics estimators agree in the low-tilt region and
> pass numerical and population-convergence checks there.  The complementary
> high-tilt region remains rare-event and boundary-term limited.  A separate
> `n=2` total-entropy endpoint estimator passes ordinary density-ratio checks
> but fails the integral-FT and exponential-weight-support gates.  Therefore
> no fluctuation-theorem symmetry is verified by the present study.

The result is not evidence that the theorem is violated.

## Next decision gate

1. The complementary `k=0.7` point and the `t=40` extension both fail the
   ancestry-support gate at `n=10`.
2. Therefore `n=20,30,40` cloning production is not launched: increasing chain
   length before fixing the rare-event support problem would not answer the FT
   question.
3. The `n=2` total-entropy control is complete and blocked by the learned
   equilibrium IFT and driven exponential-weight ESS gates.  Any continuation
   requires a better endpoint-density/rare-event method, not a larger ordinary
   sample run with the present histogram estimator.
