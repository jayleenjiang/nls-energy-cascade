# Module 06: short-chain stabilization and phase locking

## Manuscript coverage

Mechanism of stabilization in the three-mode chain.

## Numerical diagnostics

### Angular symmetry breaking

The saved NESS density is compared with its
theta1/theta3-transposed density.  The manuscript figure uses the difference
qualitatively because numerical max/mean asymmetries depend on the
high-density mask and network normalization.

### Phase-locking branch

The reduced four-dimensional system and its neural density are used to compare
observed angular concentration with the stable branch of the phase-locking
fixed-point equation.  Peak locations are diagnostics, not independent
high-precision estimators.

### Middle-mode current balance

The stationary identity balancing the two currents through the unforced middle
mode is recomputed from the saved density.  Its residual is a solver check, not
a headline result.

## Code and outputs

- Reduced-system notebook: NN notebooks/FKE_4d.ipynb.
- Five-dimensional density notebook: NN notebooks/FKE_5d_NLS.ipynb.
- Peak helper: python/analysis/find_peak.py.
- Figure: Paper/revision/symmetry_breaking.png.
- Metrics: Paper/revision/short_chain_nn_rerun_metrics.json.

## Status

The qualitative stabilization mechanism is advisor-synced.  Masked asymmetry,
peak-angle, and current-balance numbers are LOCAL_ENHANCED diagnostics and must
retain their caveats if shown.
