# Devil's Advocate review — proof note checkpoint

## Verdict: REVISE, with no fatal algebraic flaw found

### Critical issues

No critical issue was found in the coercivity, nonexplosion, heat, or
tilted-generator coefficient calculations.  The principal scientific claim is
nevertheless blocked from being called an unconditional finite-chain SCGF
theorem by the major issue below.

### Major issues

1. **The spectral conclusion still assumes the hardest existence theory.**
   The formal identity `L_k^* = Theta L_{1-k} Theta` does not by itself prove
   that a NESS exists, that the Feynman--Kac growth rate is finite, or that it
   is an isolated principal eigenvalue on a common weighted space.
   Recommendation: retain Assumption 5.2 explicitly and make the next proof
   target a model-specific Lyapunov/hypoelliptic spectral argument.

2. **The finite-time path relation needs a precise reverse ensemble.**
   With odd variables, the reverse initial density is
   `rho_t^R(w)=rho_t(Theta w)`.  A forward-only detailed relation is not
   automatic.  Recommendation: state the reverse density explicitly and keep
   the exact claim at IFT/Crooks level unless self-conjugacy is separately
   proved.  This revision has been applied to the TeX note.

3. **The right-bath gauge is an unbounded coboundary.**
   Algebraic similarity by `exp(k beta_L E)` may fail on an ill-chosen operator
   domain.  Recommendation: prove exponential-energy moments and domain
   preservation on the eventual weighted space; until then describe the gauge
   equivalence as conditional.

4. **The thermodynamic-limit theorem may be physically empty.**
   Pointwise convergence to `psi_infinity=0` preserves symmetry but does not
   establish a nontrivial fluctuation law.  Recommendation: separate existence
   of the unnormalized limit from a predeclared normalized scaling limit.

### Minor issues

- A sign typo in the enhanced manuscript's left boundary angular bath drift
  conflicts with the Cartesian SDE.  It is recorded in the module README and
  proof note; the advisor-synchronized snapshot has the correct sign.
- The symbolic script is a regression check, not evidence replacing the
  analytic derivation.
- A future manuscript insertion should cite the exact function-space theorem
  used for the Feynman--Kac principal eigenvalue rather than citing a general
  FT paper alone.

### Strongest counter-argument

> The note proves only a formal adjoint identity.  Without an
> unequal-temperature invariant law and a principal-eigenvalue theorem on the
> unbounded state space, it has not yet proved that the simulated SCGF exists,
> so it cannot promote the numerical long-chain result to an unconditional
> theorem.

This counter-argument is correct and is preserved as the central open item,
not explained away.

### Stress-test summary

| test | result |
|:---|:---|
| Remove the numerical evidence: does the generator identity remain? | yes |
| Remove the spectral assumption: does the SCGF theorem remain? | no |
| Replace finite `n` by `n -> infinity`: does the proof remain? | no |
| Could the thermodynamic limit be trivially zero? | yes |
| Is action current used as entropy production? | no |
