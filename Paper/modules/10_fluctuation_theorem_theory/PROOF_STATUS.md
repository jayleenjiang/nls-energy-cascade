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

## Two-site theorem closed in the working note

For `n=2`, both sites are directly thermostatted and the diffusion is uniformly
elliptic.  The companion proof `n2_ft_theorem.tex` establishes:

- a unique unequal-temperature NESS with a smooth strictly positive density;
- exponential NESS energy moments for every `a < 1/max(T1,T2)`;
- the finite-time forward--reverse Crooks relation and integral FT for total
  entropy;
- existence of a simple principal medium-entropy SCGF for every interior tilt
  `0 < k < 1`;
- the unconditional interior symmetry `psi_2(k)=psi_2(1-k)`.

The proof uses an exponential Foster--Lyapunov estimate, uniform ellipticity,
the Lyapunov--minorization theorem for Feynman--Kac kernels, and a common
compact-resolvent realization in `L^2(dz)`.  The model-specific algebra has an
independent symbolic check.  This is an internally closed theorem draft, not
yet an advisor- or referee-reviewed theorem.

## Not yet proved

1. Existence and uniqueness of the unequal-temperature NESS for every
   `n>=3`.
2. A model-specific Lyapunov/minorization or hypocoercive argument strong
   enough to prove that NESS result.
3. Compactness/quasi-compactness and a principal eigenvalue theorem for the
   degenerate tilted Feynman--Kac semigroup at `n>=3`.
4. A common, nonempty tilt domain with the endpoint exponential moments needed
   for the right-bath gauge.
5. Existence, uniformity, and nontrivial normalization of the
   `n -> infinity` SCGF limit.

## Allowed theorem-level wording now

For `n=2`, subject to external mathematical review of the working proof:

> The two-site continuous Cartesian NLS diffusion has a unique nonequilibrium
> steady state.  For every `0<k<1`, its medium-entropy SCGF exists and obeys
> `psi_2(k)=psi_2(1-k)`; total entropy also obeys the finite-time
> forward--reverse fluctuation relation.

For arbitrary fixed finite chain length:

> For every fixed finite chain length, the continuous Cartesian NLS diffusion
> is nonexplosive and its entropy-production tilted generator has the exact
> Gallavotti--Cohen formal-adjoint symmetry.  Consequently, any principal SCGF
> admitted by the corresponding Feynman--Kac semigroup obeys
> `psi_n(k)=psi_n(1-k)`.

The arbitrary-chain statement remains conditional at the SCGF existence step.
Neither result is an unconditional `n -> infinity` theorem.
