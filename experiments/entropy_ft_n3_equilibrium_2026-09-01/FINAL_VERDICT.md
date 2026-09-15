# Final verdict

The Cartesian bath-heat and medium-entropy accumulator **passes the frozen
equal-temperature known-answer audit** for `n=3`, `dt=0.0005`, burn-in `500`,
and block duration `20` at both `T=6` and `T=10`.

- Every requested mean (`Q_L/t`, `Q_R/t`, `J_E`, and `Sigma_m/t`) has a
  whole-stream bootstrap 95% confidence interval containing zero.
- `mean(Sigma_m/t)` is `-1.854788927375325e-06` at `T=6`, or
  `-1.8003464310168162` stream standard errors from zero.  Its 95% CI is
  `[-3.901957679884457e-06, 1.3515396393200435e-07]`.
- `mean(Sigma_m/t)` is `-1.2079929910926944e-07` at `T=10`, or
  `-0.13445406138944435` stream standard errors from zero.  Its 95% CI is
  `[-1.928062941512113e-06, 1.6848914975604784e-06]`.
- Positive/negative counts are `500097/499967` at `T=6` and
  `499962/500102` at `T=10`.
- The frozen symmetry slopes are `-1.180344402677124e-04` with 95% CI
  `[-2.739557852662919e-03, 2.499806614471104e-03]` at `T=6`, and
  `-1.257201309898014e-04` with 95% CI
  `[-2.985093886519155e-03, 2.759714954451173e-03]` at `T=10`.
  Both contain the exact equilibrium value zero.

At source level, the bath update is
`-gamma*grad(E)*dt + sqrt(2*gamma*T)*dW`; hence the continuous bath generator
annihilates `exp(-E/T)` with the same input `T`, not `2T` or `T/2`.

This audit rules out the proposed overall sign/factor-two failure at the tested
parameters and validates the normalization of the driven bath-heat and medium-
entropy observables.  It does not by itself establish a driven fluctuation
theorem, supply the system-entropy boundary term, or cure unresolved rare-tail
sampling.
