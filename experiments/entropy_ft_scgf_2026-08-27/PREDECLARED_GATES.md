# Predeclared gates

## Direct SCGF pilot

A point is marked directly reliable only if all four conditions hold:

- sample-weight ESS at least 1000;
- independent-stream ESS at least 32 of 128;
- largest individual sample contributes at most 1% of total weight;
- largest independent stream contributes at most 10% of total weight.

The stream bootstrap is paired across `k` and `1-k`.  Symmetry is judged only
where both members of a pair pass.  `k=1/2` is algebraically self-symmetric and
is not evidence for the fluctuation relation.

## Cloning gate

Implementation begins after the direct pilot identifies unsupported `k`
values.  Before any physical claim, the sampler must pass:

1. `k=0` recovery of unbiased stationary means;
2. agreement with direct SCGF where direct ESS is adequate;
3. convergence under clone-population increase;
4. convergence under observation-time increase;
5. timestep-halving consistency;
6. selection-interval convergence;
7. independent-run uncertainty and ancestry/lineage diagnostics;
8. equal-temperature and bath-swap controls;
9. no midpoint failures and the existing heat/first-law gates.

Failure of a gate leaves the result as an algorithmic pilot.  Agreement of
`psi(k)` and `psi(1-k)` would be reported as numerical consistency with the
asymptotic fluctuation relation, not as a mathematical proof.

The following quantitative rules are fixed before inspecting the `N_c=1024`
population-convergence outputs:

- use at least four independent seeds for every accepted production point;
- require minimum selection-weight ESS at least `0.1 N_c`, at least 32
  surviving initial roots, and root-count ESS at least 16;
- where direct sampling is reliable, require the cloning mean to be within
  three independent-run standard errors of the direct estimate;
- call population size converged only when consecutive populations differ by
  at most two combined standard errors and by at most `0.02` in absolute SCGF;
- apply the same two-standard-error rule, with an absolute tolerance `0.01`,
  to timestep and selection-interval comparisons;
- require zero midpoint failures in every accepted run.

These are acceptance rules, not error bars.  Passing them does not remove the
finite-time system-entropy boundary term.  A Gallavotti--Cohen comparison is
attempted only when both `k` and `1-k` separately pass all support and
convergence gates.

## Small-chain total-entropy control

The `n=2` endpoint-density study is separate from the medium-entropy cloning
study.  Its smoothing scale is selected using equal-temperature data only,
where the exact stationary density ratio is
`log rho(X_0)/rho(X_t) = Delta E/T`.  The driven data are not inspected during
this selection.

Before interpreting a learned NESS endpoint term, require all of:

- cross-fitted equal-temperature endpoint support at least 95%;
- equal-temperature learned-versus-exact log-density-ratio slope in
  `[0.9,1.1]`, correlation at least 0.9, and RMSE at most 0.5;
- the stream-bootstrap 95% interval for
  `log mean exp(-Delta s_total)` includes zero in the learned equilibrium
  control;
- exponential-weight ESS at least 1000 and at least 32 independent streams;
- maximum start/end stationarity discrepancy at most three stream-level
  standard errors across energy, both log actions, and sine/cosine phase;
- zero midpoint failures and the existing first-law residual gate.

Failure leaves the NESS total-entropy calculation exploratory.  Even after a
pass, histogram symmetry is reported only in raw two-sided bins with adequate
counts; the check is not extrapolated from `n=2` to longer chains.
