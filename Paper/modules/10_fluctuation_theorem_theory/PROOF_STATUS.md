# Proof status: boundary-driven resonant NLS fluctuation theorem

Date: 2026-08-30

## Unconditional finite-chain results established in the proof note

- The physical energy is coercive:
  `E >= (1/4) sum_j |c_j|^4 >= |z|^4/(4n)`.
- For every fixed finite `n>=2`, the Cartesian SDE has a unique global strong
  solution and is nonexplosive.
- At equal endpoint temperatures the normalized Gibbs density
  `rho_T proportional exp(-E/T)` is invariant (uniqueness is not asserted).
- The Stratonovich bath heat and its Ito form agree with the production
  accumulator:
  `dQ_r = [-gamma |F_r|^2 + gamma T_r Delta_r E]dt
          + sqrt(2 gamma T_r) F_r.dW_r`, with `Delta_r E=4M`.
- On `C_c^infinity(R^{2n})`, the medium-entropy tilted generator satisfies the
  exact formal-adjoint identity
  `L_{n,k}^* = Theta L_{n,1-k} Theta`.

## Conditional theorems established

- If the long-time Feynman--Kac growth rate is the spectral bound/principal
  eigenvalue of the closed tilted generator on a common function space, then
  `psi_n(k)=psi_n(1-k)` for that finite `n`.
- The right-bath gauge used in the cloning experiment is spectrally equivalent
  to medium entropy whenever the energy-coboundary conjugation preserves the
  chosen operator domains.
- If `psi_n(k)` exists, has the finite-chain symmetry, and converges pointwise
  to `psi_infinity(k)` on a symmetric tilt interval, then the limiting SCGF
  inherits the same symmetry.

## Not yet proved

1. Existence and uniqueness of the unequal-temperature NESS for every `n`.
2. A model-specific Lyapunov/minorization or hypocoercive argument strong
   enough to prove that NESS result.
3. Compactness/quasi-compactness and a principal eigenvalue theorem for the
   tilted Feynman--Kac semigroup on the required unbounded state space.
4. A common, nonempty tilt domain with the endpoint exponential moments needed
   for the right-bath gauge.
5. Existence, uniformity, and nontrivial normalization of the
   `n -> infinity` SCGF limit.

## Allowed theorem-level wording now

> For every fixed finite chain length, the continuous Cartesian NLS diffusion
> is nonexplosive and its entropy-production tilted generator has the exact
> Gallavotti--Cohen formal-adjoint symmetry.  Consequently, any principal SCGF
> admitted by the corresponding Feynman--Kac semigroup obeys
> `psi_n(k)=psi_n(1-k)`.

This is stronger than a numerical fit but remains conditional at the SCGF
existence step.  It is not yet an unconditional `n -> infinity` theorem.
