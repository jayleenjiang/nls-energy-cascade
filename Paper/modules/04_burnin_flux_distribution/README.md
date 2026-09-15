# Module 04: burn-in and finite-time current distribution

## Manuscript coverage

Latest mentor questions on relaxation, burn-in selection, and the dependence of
log P[J_tau>A] on averaging time tau and threshold A.

## Experiment 1: no-burn-in relaxation

- Chains: n=10,20,40,80.
- Checkpoints: T=50,100,...,500.
- Records: terminal mean action profile, cumulative profile, cumulative current,
  and last-interval current.
- Purpose: choose burn-in from observable relaxation rather than from a fixed
  arbitrary number.

## Experiment 2: finite-time current tails

- Chains: n=10,20,30,40.
- Averaging times: tau=20,40,...,200.
- Threshold grid: A=0.01,0.02,... in the analysis.
- Common pilot burn-in: 500.
- Target sample size: 100000 trajectories per n.
- Simulator: flux/NLS_flux_relaxation_tau.cpp.
- Analysis: flux/analyze_burnin_ld.py.

The analysis sorts the trajectory currents, computes empirical survival
probabilities, plots log P[J_tau>A] versus A, and studies its tau dependence.
An approximately straight upper-middle segment is descriptive evidence only;
large-deviation scaling requires testing whether -log P divided by tau
collapses as tau grows.

## Run state

The pilot with 1024 trajectories is complete.  The full 100k run writes to
Paper/revision/experiments/flux_validation/burnin_ld_full100k_2026-08-25 and
was active when this module was created.  Do not move or rename that directory.
After all four n values finish, run the registered analysis script in the same
directory and update this README from RUNNING to ANALYZED.
