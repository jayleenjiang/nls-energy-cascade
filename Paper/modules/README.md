# NLS paper research modules

Updated: 2026-08-25

This directory is the authoritative navigation layer for the paper.  It does
not replace or move the historical working directories.  Every module records:

1. the manuscript section and scientific question;
2. the canonical simulation or analysis code;
3. raw and processed data;
4. manuscript-ready figures and numerical reports; and
5. whether the material has been synchronized with the advisor.

The scientific observable called flux in older filenames is the conserved
action current unless a module explicitly states otherwise.  It is not
automatically the Hamiltonian heat current.

## Module map

| ID | Manuscript content | Current status |
|---|---|---|
| 00 | Manuscript versions and build packages | Two paper lines kept separate |
| 01 | Introduction, derivation, cascade geometry, thermodynamic setup | Theory/source-traced geometry |
| 02 | Long-chain NESS profiles and local thermodynamic equilibrium | Core content advisor-synced; residual mesh advisor-requested |
| 03 | Action-current conductivity and finite-size scaling | Canonical production result; local enhanced text needs advisor approval |
| 04 | Burn-in selection and finite-time flux distribution | Latest advisor-requested experiment; full 100k run active |
| 05 | Three-mode Fokker--Planck NESS density | Core short-chain content advisor-synced |
| 06 | Short-chain stabilization and phase locking | Qualitative core synced; quantitative rerun diagnostics local |
| 07 | Backward-generator relaxation and eigenfunction surrogate | Local enhanced diagnostic; not a resolved spectral gap |
| 08 | Numerical validation and robustness appendix | Validation-only extensions; do not silently migrate into advisor draft |
| 09 | Reproducibility, audits, packaging, and submission materials | Local workflow and author-confirmation layer |

## Status vocabulary

- ADVISOR_SYNCED: already present in the advisor-facing scientific baseline.
- ADVISOR_REQUESTED: explicitly requested in the latest meetings.
- LOCAL_ENHANCED: locally added analysis or manuscript text not yet approved.
- VALIDATION_ONLY: useful robustness evidence, but not part of the synchronized
  main narrative.
- RUNNING: outputs are still being written and must not be moved.
- LEGACY_REFERENCE: retained only to explain historical comparisons; not a
  source for current claims.

## How to use this directory

- Read MODULE_INDEX.tsv for a one-row-per-module overview.
- Read FILE_ORIGINS.md before deciding which duplicate-looking directory is
  canonical.
- Each module has README.md and manifest.tsv.
- Run build_module_view.sh to create local/code, local/data, local/figures, and
  local/reports symlinks.  These links are ignored by Git and never duplicate
  the underlying data.
- Run check_module_links.py after adding or moving any source artifact.

## Preservation rules

- No historical code, data, figure, or manuscript file was deleted or moved in
  this organization pass.
- The pre-organization root README and the two compared manuscripts are backed
  up under Paper/backups/2026-08-25_repo_organization.
- Paper/revision/experiments/flux_validation/burnin_ld_full100k_2026-08-25 is
  an active output directory and remains in place until the run and analysis
  finish.
- Paths containing Paper/revision_2026-06-19 in older reports are historical.
  The live local revision root is Paper/revision.

## Update protocol

When a new experiment is added:

1. put executable source in the existing canonical source tree;
2. write outputs to one dated experiment directory;
3. add that directory and its analysis script to exactly one primary module;
4. add cross-links from validation modules only when needed;
5. mark advisor status explicitly; and
6. run the link checker before updating manuscript claims.
