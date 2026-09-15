# Final verdict

This analysis used only the saved stationary correlation arrays. No simulation
was started.

The requested five-observable even-candidate group has no common plateau in
either bath condition. In the driven case only `I2` and `I1+I3` resolve
individual plateaus, fewer than the frozen minimum of three. Their rates span
`[-1.109,-0.993]` (11.1% relative spread). At equal temperature,
`cos(theta1)`, `cos(theta3)`, and `I2` resolve separately, but their rates span
`[-1.188,-0.907]` (24.8%) and their 95% intervals do not intersect.

On the predeclared early window `[0.05,1.00]`, `I1-I3` has a resolved damped
mode in both cases:

- driven: `lambda_R=-4.911 [-5.146,-4.694]`,
  `lambda_I=5.356 [5.205,5.507]`;
- equal temperature: `lambda_R=-4.285 [-4.427,-4.147]`,
  `lambda_I=5.351 [5.200,5.499]`.

`sin(theta1)` gives a distinct damped fit. It resolves at equal temperature,
`lambda_R=-5.006 [-5.069,-4.945]`, `lambda_I=3.896 [3.714,4.066]`, but the
driven estimate `lambda_I=3.306 [3.112,3.495]` misses the frozen half-cycle
threshold `3.30694` by `9.1e-4` and is therefore marked unresolved.

The two sign-changing observables have disjoint confidence intervals for both
the real and imaginary parts in both bath conditions. They do not identify one
common complex-conjugate eigenpair. Their fitted zero crossings agree closely
with the observed crossings, so the damped description is useful, but a sign
change plus preference over a single real exponential is not by itself a proof
of a complex generator eigenvalue.

Finally, `sin(theta1)` alone is not exchange odd: endpoint exchange maps it to
`sin(theta3)`, which was not saved. Exact sector resolution requires the
left/right sums and differences, especially `sin(theta1) +/- sin(theta3)`.
