# Module 10: fluctuation-theorem proof program

## Question

For the projection-free Cartesian boundary-driven resonant NLS chain, which
parts of the finite-chain and thermodynamic-limit fluctuation theorem can be
proved from the continuous-time model, independently of the numerical
integrator?

## Frozen continuous model

The proof uses the physical energy `E=H/2` and the Cartesian SDE implemented in
`flux/NLS_entropy_ft.cpp` and `flux/NLS_entropy_cloning.cpp`.  It does not use
the projected action--angle Euler scheme.

## Files

- `finite_n_ft_proof.tex`: theorem/lemma/proof note with an explicit proof
  ledger.
- `PROOF_STATUS.md`: concise statement of what is proved and what remains an
  analytic assumption.
- `verify_algebra.py`: independent symbolic checks of the energy gradient,
  time-reversal invariance, boundary Laplacian, and tilted-generator
  coefficient identity for a representative chain.

## Current status

The following are proved for every fixed finite `n>=2`:

1. coercivity of the physical energy;
2. global nonexplosion of the Cartesian diffusion;
3. invariance of the equal-temperature Gibbs law;
4. the medium-entropy tilted-generator identity
   `L_k^* = Theta L_{1-k} Theta` on compactly supported smooth test functions;
5. inheritance of the symmetry by any existing finite-chain SCGF;
6. inheritance by any existing thermodynamic-limit SCGF.

The remaining hard steps are existence/uniqueness of the unequal-temperature
NESS, the required Feynman--Kac spectral theorem on an unbounded state space,
and existence/nontriviality of the `n -> infinity` SCGF limit.  The note does
not label those steps as proved.

## Model-consistency finding

The Cartesian bath implies the left-end action--angle drift

```text
b_phi_1 = -2 gamma I_2 sin(theta_2),
```

and the right-end drift

```text
b_phi_n = +2 gamma I_{n-1} sin(theta_n).
```

`Paper/revision/advisor_safe_migration_2026-08-26/draft_advisor_safe.tex`
has these signs.  The enhanced `Paper/revision/draft.tex` currently displays a
plus sign at the left end and must not be used as the frozen theorem statement
until that typo is corrected with a manuscript backup.
