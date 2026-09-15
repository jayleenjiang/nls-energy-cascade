# Final verdict: Section 4.2 stabilization rebuild

Analysis completed without a new simulation or density retraining.

## Outcome

The frozen all-diagnostic verdict is **FAIL**.  This is not a failure of every
physical conclusion; it means that the validated density does not pass every
predeclared Section 4.2 cross-check.

- Endpoint exchange: the driven angular-TV point estimate is 4.002 times the
  measured equilibrium point floor and all three density seeds agree tightly.
  The density value 0.11489 also agrees with the direct 0.11587.  The stricter
  frozen lower-bound / three-times-upper-floor gate is 2.435 rather than 3 and
  therefore fails.
- Phase locking: 27 of 30 new-density versus held-out interval comparisons
  pass.  Three fail: low-energy theta3 circular mean, middle-energy theta1
  circular mean, and high-energy theta1 locking probability.
- Current balance: the complete held-out driven trajectories satisfy the exact
  global identity, with J12+J23 = 1.7046e-5 and 95% CI
  [-2.7013e-4, 3.1681e-4].  On the Section 4.1 support mask, however, the new
  density and trajectory estimates differ by -2.306 bootstrap standard errors,
  so the density reconstruction gate fails.

## Claim boundary

The updated evidence supports the qualitative stabilization picture on the
approximately 99.4% sampled-support region: reproducible exchange asymmetry,
energy-dependent phase locking, and global direct current balance.  It does
not support a claim that every Section 4.2 statistic is quantitatively
reproduced by the density.  The three phase failures and masked-current
discrepancy must accompany any paper-level statement.

The legacy density is retained in every raw comparison table.  It is generally
farther from direct trajectories and has only 3.37% effective sample fraction
in the driven importance-sampling calculation.
