# Theory and interpretation guardrails

This note fixes what the numerical experiment can and cannot claim before the
large production sample is inspected.

## 1. Three distinct random variables

The joint sampler stores three related but noninterchangeable quantities.

1. The bath heat current

   \[
   J_E(t)=\frac{Q_L-Q_R}{2t}.
   \]

2. The medium entropy production

   \[
   \Sigma_t^{\mathrm m}=-\frac{Q_L}{T_L}-\frac{Q_R}{T_R}.
   \]

3. The bulk action current

   \[
   J_M(t)=\frac1t\int_0^t j_M(s)\,ds.
   \]

Only the first two are directly thermodynamic.  The action-current symmetry is
an empirical finite-time rate diagnostic unless an additional derivation shows
that it is tightly coupled to the entropy current.

For the present physical-energy convention `E=H_code/2=H_paper`, the Gibbs
weight is `exp(-E/T)`, so `1/T` is the inverse temperature multiplying bath
heat.  The complex-Langevin stochastic-thermodynamics treatment likewise
separates entropy, heat, particle/action, and energy currents; see Borlenghi et
al., [Entropy production for complex Langevin
equations](https://arxiv.org/abs/1704.01566).

## 2. Exact finite-time balance identity

Let

\[
D_E(t)=\frac{\Delta E}{t}=\frac{Q_L+Q_R}{t}+
O(\text{integrator error}).
\]

Then the stored observables obey

\[
\frac{\Sigma_t^{\mathrm m}}{t}
=\left(\frac1{T_R}-\frac1{T_L}\right)J_E(t)
-\frac12\left(\frac1{T_L}+\frac1{T_R}\right)D_E(t)
+O(\text{integrator error}).
\]

For `T_L=10` and `T_R=2`, the heat-current affinity is `0.4` and the
finite-time energy-boundary coefficient is `0.3`.  In the stationary mean,
`<D_E>=0`, giving `<Sigma_m/t> = 0.4 <J_E>`, which is satisfied by the pilot
for every chain length.

## 3. Which fluctuation theorem is exact

The finite-time integral/detailed fluctuation theorem applies to total
trajectory entropy,

\[
\Delta s_{\rm tot}=\Sigma_t^{\mathrm m}+\Delta s_{\rm sys},
\qquad
\Delta s_{\rm sys}=-\log\rho(X_t)+\log\rho(X_0),
\]

with the appropriate initial/final ensemble and time-reversal operation.  The
system-entropy term is not available from bath heat alone.  This distinction is
standard in stochastic thermodynamics; see Seifert,
[Entropy production along a stochastic trajectory and an integral fluctuation
theorem](https://arxiv.org/abs/cond-mat/0503686).

The complex amplitudes also contain even and odd components under conjugation
time reversal.  The parity-aware path-probability construction is discussed by
Spinney and Ford, [Non-equilibrium thermodynamics of stochastic systems with
odd and even variables](https://arxiv.org/abs/1201.0904).  In the present
Cartesian representation, `x` is even and `y` is odd; the Hamiltonian drift is
reversible and the gradient bath drift is irreversible under this parity.

At long times, the medium entropy and heat-current symmetry may approach the
Gallavotti--Cohen form because endpoint contributions are nominally `O(1)`.
However, unbounded energy fluctuations can make heat boundary terms affect
rare tails even asymptotically.  Extended heat fluctuation relations are known
in such Langevin systems; see van Zon and Cohen, [An Extension of the
Fluctuation Theorem](https://arxiv.org/abs/cond-mat/0305147).

## 4. How the production data will be judged

The following hierarchy is predeclared.

1. **Raw symmetric-bin diagnostic.**  Fit

   \[
   R_t(a)=t^{-1}\log[p_t(a)/p_t(-a)]
   \]

   only in bins with nonzero raw counts on both sides and at least the specified
   effective count.  Plus-four smoothing is visual only.

2. **Time dependence.**  A slope close to one at a single `t` is insufficient.
   The entropy slope must be stable or move systematically toward one as `t`
   increases, with uncertainty and usable overlap reported.

3. **Heat-current cross-check.**  The corresponding raw heat-current symmetry
   has reference slope

   \[
   \Delta\beta=T_R^{-1}-T_L^{-1}=0.4.
   \]

4. **Normal tail benchmark.**  A Gaussian CDF is compared with each empirical
   tail, but its fitted slope is not an FT test.  The FT concerns a probability
   ratio; Gaussianity is neither necessary nor sufficient.

5. **Action-current result.**  Report the empirical action-current symmetry
   slope without labeling it entropy production.  Correlation with `J_E` is
   accompanied by the regression residual variance; correlation alone does
   not prove tight coupling.

6. **Numerical gates.**  No physical conclusion is accepted unless the run has
   zero midpoint failures, no nonfinite values, stationary means, converged
   timestep subsets, and a small `Q_L+Q_R-DeltaE` residual.

## 5. Allowed conclusions

- If the raw entropy/heat symmetry slopes agree with their references over a
  stable range of `t` and `a`, the result is **consistent with** the asymptotic
  fluctuation relation in that tested window.
- If they do not, the correct conclusion is **not verified in the sampled
  window**, not that the theorem is universally false.
- If negative events disappear at larger `t`, direct sampling is
  resolution-limited.  A later tilted-generator/cloning or other importance
  sampling study would then be a separate experiment, not an extrapolation of
  zero-count plus-four bins.
