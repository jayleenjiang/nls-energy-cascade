# Module 03: action-current scaling

## Manuscript coverage

Action-current conductivity and finite-size transport scaling.

## Canonical production experiment

- Scientific question: how does the stationary conserved action current depend
  on chain length?
- Source: flux/NLS_flux_canonical.cpp.
- Parameters: (T1,Tn,gamma)=(10,2,0.1), fixed-step Euler--Maruyama,
  dt=5e-4, measurement window 200.
- Chains: n=10,20,30,40.
- Independent units: 1024 trajectories per n.
- Burn-ins: 1000,1280,2880,5120.
- Primary result directory:
  Paper/revision/experiments/flux_validation/production_dt5e-4.
- Reported fit: E[J(n)] = 28.7457 n^-1.85008, log-fit R2=0.998013,
  trajectory-bootstrap 95% CI [-1.87034,-1.83049].

This is a finite-size law for the action current.  It is not by itself an
asymptotic theorem and should not be called Hamiltonian heat conductivity.

## SIMD comparison

The dated fixed-SIMD dataset and NLS_flux_SIMD.cpp are retained to compare two
implementations of the same SDE.  They are not the primary source of the paper
exponent.  Older flux_V1/flux_V2 results and pre-fix datasets are
LEGACY_REFERENCE only.

## Status

The canonical numbers are source-traced production outputs, but the complete
conductivity section was added in the local enhanced draft after the
advisor-facing draft still contained a TODO.  Present the four production
means and fit to the advisor before treating the section as synchronized.
