# Final verdict

The frozen five-dimensional KDE route does **not** reliably resolve the
finite-time total-entropy fluctuation theorem at `n=3`, `t=20`, and
`(T_L,T_R)=(10,2)`.

This is a negative numerical-method result, not evidence that the physical
fluctuation theorem is false.

- The equilibrium Gibbs known-answer KDE gate fails at both `T=6` and `T=10`.
- The two independently trained driven KDEs disagree beyond every frozen
  accuracy threshold.
- The frozen binned KDE has zero support for 42 of 1,000,064 driven endpoint
  pairs.  Per the no-extrapolation rule, total-entropy sign counts, DFT slope,
  and IFT were not computed.
- The medium-entropy result remains valid: 41 negative blocks out of
  1,000,064, but zero symmetric bin pairs contain 20 observations on both
  sides.  Its exponential average has ESS 2.288 and is unresolved.

Therefore no paper-level claim that the exact finite-time total-entropy FT has
been extended from `n=2` to `n=3` is supported by this KDE experiment.
