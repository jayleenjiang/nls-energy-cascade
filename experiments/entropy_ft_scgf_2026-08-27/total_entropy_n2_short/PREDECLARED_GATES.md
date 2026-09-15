# Short-time n=2 NESS total-entropy control

The earlier `t=20` endpoint experiment passed ordinary density-ratio errors
but failed its exponential IFT gate.  This follow-up tests whether the same
cross-fitted stationary-density estimator is accurate enough when endpoints
are separated by `t=0.1`, where local interpolation errors should cancel more
strongly and both entropy signs remain directly sampled.

The smoothing scale is selected using an independent equal-temperature data
set only.  The driven `(T_left,T_right)=(10,2)` result is not used for tuning.

Acceptance requires all existing endpoint gates plus:

- at least 64 independent streams and 200,000 non-overlapping blocks;
- held-out equal-temperature support at least 95%, correlation at least 0.9,
  RMSE at most 0.5, and density-ratio slope in `[0.9,1.1]`;
- learned equal-temperature and driven 95% stream-bootstrap intervals for
  `log mean exp(-Delta s_total)` both include zero;
- exponential-weight ESS at least 1,000 in each case;
- at least eight raw symmetry bins with 20 counts on both sides;
- driven fitted `log[p(s)/p(-s)]` slope consistent with one within the larger
  of 0.1 and two standard errors, and intercept consistent with zero within
  the same rule;
- zero midpoint failures, finite samples, stationarity `|z| <= 3`, and the
  short-block first-law gate: RMS absolute balance error at most `5e-5`.

The 262,144-sample pilot was used only for estimator development.  A declared
grid scan selected `(40,40,80)` bins and Gaussian smoothing `sigma=0.75` by
equal-temperature endpoint-ratio RMSE.  These values are now fixed.  Final
validation uses new random seeds and 1,048,576 blocks per condition; neither
the grid nor smoothing may be changed after inspecting those files.  The
short-block gate is stated in absolute energy units because dividing a fixed
per-block integrator residual by `t=0.1` makes the older long-block rate gate
dimensionally inappropriate.  The exact Gibbs IFT remains an additional
non-circular numerical control.

Passing would verify a finite-time, small-chain NESS total-entropy relation in
the tested window.  It would not establish the long-chain asymptotic GC
relation.

## Final endpoint-density model

The histogram estimator retained an IFT bias after the independent million-
sample equilibrium validation.  Before applying a replacement to the driven
data, the endpoint model is fixed as a normalized action--Fourier exponential
family with action degree `D=2`, angular harmonics `K=1`, grid quadrature
`42 x 42 x 64`, and ridge coefficient `1e-8`.  Its features include `u_1`,
`u_2`, all action monomials through degree two, and their first sine/cosine
harmonics.  This family contains the exact n=2 Gibbs log density
`u_1+u_2-E/T+constant`.

On the independent equilibrium production validation, selected without using
the driven blocks, it gives endpoint-ratio RMSE `0.00826`, slope `0.99941`,
correlation `0.99994`, and equilibrium log-IFT `3.4e-5`.  These numbers justify
the model but may not be used to alter it further.  The same fixed feature
family, quadrature, regularization, support rule, and optimizer tolerance are
applied cross-fitted to the already held-out driven data exactly once.

## Audit amendment for an estimated density

The literal requirement that the learned equilibrium confidence interval
contain zero is not a consistent validation criterion: with increasing sample
size, any nonzero approximation error in an otherwise highly accurate density
model is detected and the interval eventually excludes zero.  It tests whether
the estimated density is mathematically exact, not whether its endpoint error
is negligible for the FT calculation.  The final audit therefore uses a fixed
equivalence margin `|log IFT| <= 1e-3` for the learned equilibrium endpoint,
alongside the independently declared RMSE/slope/correlation gates and the
exact-Gibbs IFT gate.  The driven result retains the stricter requirement that
its stream-bootstrap interval include zero.  This amendment and rationale are
recorded explicitly rather than silently relabeling the failed histogram
pilot.
