# Module 02: long-chain NESS profiles and LTE

## Manuscript coverage

Long-chain nonequilibrium steady state, action profiles, local thermodynamic
equilibrium, and residual diagnostics.

## Experiments

### A. Steady action profiles

- Chains: n=25, 50, 100 with endpoint temperatures 10 and 2.
- Purpose: establish a stabilized monotone action profile and define local
  kinetic-temperature estimates.
- Primary simulator family: cpp/simulation/lte_histogram_simd.cpp and
  lte_simd_dt25.cpp.
- Profile inputs: experiments/lte/test_profile.txt,
  experiments/lte/simd_n50_profile.txt, and
  experiments/lte/n100_dt25_profile.txt.
- Figure: Paper/revision/action_profiles.pdf.

### B. Local-equilibrium pair histograms

- Observable: three-dimensional pair histogram of
  (I_j, I_{j+1}, theta_j).
- Chains/sites: n=25, 50, 100 at normalized positions near 0.24, 0.48, 0.72.
- Purpose: compare the NESS pair marginal q with simulated equilibrium
  marginals p_T and quantify log q = x log p_T + c.
- Raw data: lte/n25 data, lte/n50 data, lte/n100 data, and temperature controls.
- Analysis: python/analysis/hist_analysis.py,
  python/analysis/lte_full_analysis.py, and source-trace scripts.

### C. Advisor-requested residual mesh

- Purpose: inspect where the rescaled-equilibrium approximation fails in the
  full three-dimensional density, rather than relying on a one-dimensional
  residual plot.
- MATLAB source preserved here: code/matlab/compare_residual.m,
  residual.m, and fit_new.m.
- Primary n=15 inputs: lte/n15 data.
- Paper-ready outputs: compare_residual_mesh.pdf and
  compare_residual_matlab.pdf.

### D. Conditioned LTE check

- Purpose: test conditional pair-density structure at fixed total action.
- Code/data: lte/cond_hist and cpp/simulation/lte_simd_conditioned.cpp.
- Role: diagnostic support, not a separate headline claim.

## Status

The profile and core LTE comparison are advisor-synced.  The n=15 3D residual
mesh was explicitly requested later.  Residual-decomposition tables and
source-trace metrics are local enhancements and should be explained before
being inserted wholesale.
