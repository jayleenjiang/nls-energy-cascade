# Predeclared gates for the discrete path-ratio control

The 20,000-trajectory runs in `pilot/` are exploratory and select only a
short-time region where direct exponential averaging is feasible.  Production
gates are fixed here before increasing the sample size or scanning timestep.

For each forward/reverse ensemble require:

- zero midpoint failures and no nonfinite trajectories;
- the 128-group bootstrap 95% interval for
  `log mean exp(-Sigma_total_discrete)` contains zero;
- exponential-weight ESS at least 1,000;
- maximum single-trajectory weight fraction at most 1%.

For the bidirectional Crooks histogram require:

- identical symmetric bins for forward `Sigma_F` and reverse `-Sigma_R`;
- no pseudocounts in the accepted fit;
- at least 50 raw samples on both sides of each accepted bin;
- at least eight accepted bins;
- fitted slope consistent with one within the larger of 0.1 and two standard
  errors;
- fitted intercept consistent with zero within the larger of 0.1 and two
  standard errors.

The heat comparison is a convergence diagnostic, not part of the exact
discrete identity.  The mean and RMS difference between the transition-kernel
medium entropy and `-Q_left/T_left-Q_right/T_right` must decrease under
timestep halving before the discrete result is connected to the continuous
physical entropy.

Passing these gates verifies a finite-time transient fluctuation relation for
the explicitly defined discrete forward/reverse path measures.  It does not by
itself establish the long-time same-process Gallavotti--Cohen symmetry in the
NESS.
