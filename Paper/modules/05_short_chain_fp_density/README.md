# Module 05: three-mode Fokker--Planck NESS density

## Manuscript coverage

Detailed property of the short chain; numerical solution of the
nonequilibrium steady-state Fokker--Planck equation.

## Experiment

- System: n=3 reduced to
  (I1,I2,I3,theta1,theta3).
- Purpose: reconstruct the stationary five-dimensional density with a
  data-assisted neural Fokker--Planck solver.
- Monte Carlo / NC-KDE source family: cpp/fp5d.
- Neural solver: NN notebooks/FKE_5d_NLS.ipynb.
- Archived training/evaluation workspace: KDE/4:15_NN.
- Equilibrium control: compare network slices to the exact Gibbs density.
- Nonequilibrium output: angular density slices showing tilt and loss of
  exchange symmetry under T1 != T3.

## Primary outputs

- Paper/revision/eq_validation.png.
- Paper/revision/neq_density.png.
- Paper/revision/short_chain_nn_rerun_metrics.json.

## Status

The equilibrium/NESS density narrative is part of the advisor-synced
short-chain core.  The exact percentages in the rerun metrics are local
source-trace diagnostics; use them only with their slice and mask definitions.
