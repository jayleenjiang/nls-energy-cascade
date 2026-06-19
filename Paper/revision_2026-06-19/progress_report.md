# Progress report — 2026-06-19

## Completed in this revision pass

1. Preserved the original extracted manuscript/code and created a revision
   working copy under `Paper/revision_2026-06-19/`.
2. Audited the draft and identified P0 model/data consistency issues.
3. Replaced the invalid June 18 flux workflow with a canonical
   Gibbs-preserving action-current simulator:
   `flux/NLS_flux_canonical.cpp`.
4. Added reproducible analysis scripts:
   - `flux/analyze_canonical_flux.py`
   - `flux/analyze_current_windows.py`
   - `flux/gibbs_mcmc_reference.py`
   - `flux/compare_gibbs_sde.py`
5. Ran validation gates:
   - warning-clean optimized C++ build;
   - sanitizer smoke test;
   - deterministic one-thread/two-thread reproducibility;
   - equal-temperature zero-current check;
   - hot/cold reversal check;
   - independent Gibbs MCMC profile check.
6. Ran production current-scaling experiments for `n=10,20,30,40` with
   corrected boundary noise and correct standard errors.
7. Updated `draft.tex` to:
   - correct the full-coordinate diffusion coefficient;
   - separate full phase variables from doubled bond-angle variables;
   - fix the Fourier-law sign convention;
   - use action-current terminology instead of heat-current terminology;
   - remove unsupported entropy/fluctuation-theorem claims;
   - replace the invalid `n^-1.75` claim with the corrected `n^-1.850`
     production result;
   - fix the stable phase-locking branch;
   - weaken unsupported “real spectrum/overdamped” claims;
   - insert the corrected scaling figure.
8. Restored local image resources needed by the short-chain figures and compiled
   the revision cleanly.
9. Added submission-critical manuscript sections:
   - English abstract;
   - conclusion and limitations;
   - data/code availability;
   - declarations for ethics, author contributions, competing interests,
     funding status, and AI-assisted preparation.
10. Replaced the remaining long-chain figure placeholders with reproducible
    generated figures:
    - `action_profiles.pdf/png`;
    - `cascade_embedding.pdf/png`;
    - `lte_residual_midchain.pdf/png`.
    The generation script is
    `Paper/revision_2026-06-19/scripts/generate_manuscript_figures.py`, and
    source metrics are recorded in
    `Paper/revision_2026-06-19/manuscript_figure_metrics.json`.
11. Replaced the phase-locking “derivation to be supplied” note with the
    reduced fixed-point quadratic
    `4(1+gamma^2) sin(theta)^2 + 4 gamma sin(theta) - 3 = 0` and the stable
    branch selection.
12. Recompiled the manuscript with TeX Live. The current PDF has 17 pages and
    the LaTeX log has no unresolved citations/references and no overfull or
    underfull box warnings.
13. Resolved two remaining rigor issues in the short-chain section:
    - removed the unsupported quantitative “about 15%” symmetry-breaking claim
      and made the figure caption explicitly qualitative;
    - added a source-traced eigen relaxation-rate sensitivity analysis in
      `scripts/analyze_eigen_fit_windows.py`, with results in
      `eigen_fit_sensitivity.json`, and replaced `eigenvalue_scatter.png`
      with a reproducible relaxation-rate diagnostic figure. The historical
      `lambda_R=-0.934` value is now described as an observable-dependent
      diagnostic rate, not a high-accuracy spectral-gap estimate.
14. Converted the reference list from an inline `thebibliography` block to
    `references.bib`, verified all eight references against publisher/arXiv/DOI
    metadata, updated `DobsonLiZhai` to its published CMS version, and recorded
    the audit in `integrity_audit_2026-06-19.md`.

## Key validated numerical result

For the corrected Gibbs-preserving SDE at `T1=10`, `Tn=2`, `gamma=0.1`,
`dt=5e-4`, burn-ins `(1000,1280,2880,5120)`, and `1024` trajectories per chain
length:

| n | mean action current | SE |
|---:|---:|---:|
| 10 | 0.3925219606 | 0.0018902080 |
| 20 | 0.1191693526 | 0.0009195305 |
| 30 | 0.0545731139 | 0.0006191849 |
| 40 | 0.0297475540 | 0.0004827205 |

Power-law fit:

- `E[J(n)] = 28.7457 n^-1.85008`
- `R^2 = 0.998013`
- bootstrap 95% exponent CI: `[-1.87034, -1.83049]`

## Remaining blockers before journal submission

- Funding statement still needs author confirmation. The current manuscript
  explicitly says funding information was not supplied in the available
  materials.
- Author-contribution and competing-interest declarations should be confirmed
  by the authors before submission.
- The complete LTE table (`tab:lte`) needs a tracked regeneration/export file;
  the current mid-chain residual figure is source-traced, but all table entries
  are not yet captured in a unified metrics artifact.
- Short-chain neural-network figure metrics (`eq_validation`, `neq_density`,
  `symmetry_breaking`, `Q1_slices`) should be exported from notebooks/scripts to
  a tracked metrics file before final submission.
- Perform the final 100% citation/data/claim audit after the LTE and
  short-chain NN metrics are exported.
- Optional but recommended: add an appendix or supplementary note documenting
  LTE histogram generation and current-scaling validation in a compact
  reproducibility table.
