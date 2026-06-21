# NLS numerical paper: material inventory

Date: 2026-06-19

## Canonical manuscript inputs

- Framework supplied by the user: `Paper/NLS_numerical_study_source_2026-06-19/paper_draft_1.tex`
- Working manuscript supplied by the user: `Paper/NLS_numerical_study_source_2026-06-19/draft.tex`
- Revision copy: `Paper/revision_2026-06-19/draft.tex`
- Original backups: `Paper/backups/2026-06-19/`
- Source archive retained unchanged:
  `/Users/jayleenjiang/Downloads/NLS_numerical_study.zip`
- Archive SHA-256:
  `bd01b3b521c67d78308d1158f09b2d526190930b8fd52ef9d6869a546cb2c41a`

The archive contains 45 files (about 31 MB uncompressed): 12 TeX files,
figures, four flux datasets, LTE code and reports, and backward Monte Carlo
training data. The framework filename in the archive is `paper_draft_1.tex`,
although it was described as `paper_draft.tex`.

## Research objective reconstructed from the materials

The paper studies a stochastic energy-cascade chain derived from the resonant
cubic NLS toy model. Its intended contributions are:

1. formulate and numerically validate a Gibbs-preserving two-bath version of
   the chain;
2. characterize long-chain nonequilibrium steady states, especially action
   profiles and approximate local equilibrium;
3. determine the chain-length scaling and finite-time distribution of the
   action current;
4. solve the three-mode stationary Fokker--Planck equation with a
   data-assisted neural network and validate it at equilibrium;
5. examine phase-dependent stabilization in the short chain; and
6. estimate slow relaxation modes/eigenfunctions of the backward generator.

## Local primary references read

- `Paper/CKSTT.pdf`: derivation and deterministic energy-transfer toy model.
- `Paper/Non-equilibrium steady state for a three-mode energy cascade model.pdf`:
  rigorous three-mode NESS result for a different, proof-oriented thermostat.
- `Paper/Summary_NLS.pdf`: derivation of the natural Gibbs-preserving bath in
  complex and action-angle coordinates.
- `Paper/zhai22a-3.pdf`: mesh-free, data-assisted neural Fokker--Planck solver.
- `Paper/1390904.pdf`: fast multivariate binned kernel estimation.
- `Paper/7cases_Notes.pdf`: historical boundary-condition experiments.

## Code/data provenance map

| Manuscript result | Principal local source | Status |
|---|---|---|
| Long-chain LTE histograms | `cpp/simulation/lte_histogram_simd.cpp` and `lte/` | Canonical corrected bath; data present |
| Conditioned LTE test | `lte/cond_hist/lte_simd_conditioned.cpp` | Corrected bath; data present |
| 5D equilibrium/NESS density | `cpp/fp5d/`, `NN notebooks/FKE_5d_NLS.ipynb`, `KDE/` | Corrected bath; requires metric audit |
| Stabilization/phase locking | `Paper/Summary_NLS.pdf`, 4D notebook/code, manuscript figures | Formula/branch and claim audit required |
| Generator relaxation/eigenfunction | `cpp/backward/NLS_backward.cpp`, `python/data_gen/get_train.py`, `NN notebooks/FKE_eigen.ipynb` | Internally inconsistent estimates |
| Flux distribution, 2026-06-18 | `flux/NLS_flux_SIMD.cpp`, `flux_n*.txt` | Invalid for the stated model: pre-fix noise |
| Older flux scaling | `cpp/simulation/flux_V1.cpp`, `flux_V2.cpp`, `Energy Cascade/Flux/` | `flux_V2` has a documented phase-drift sign bug |

The archive copies of `NLS_backward.cpp` and `lte_histogram_simd.cpp` are
byte-identical to the corresponding repository copies (verified by SHA-256).

## Build and visual baseline

- `draft.tex` compiles with TeX Live to a 13-page PDF.
- All cross-references resolve.
- There is one overfull paragraph around the LTE effective-temperature text.
- Four long-chain figures are placeholders.
- The compiled manuscript contains no abstract, no conclusion, no
  reproducibility/data-availability statement, and no appendix documenting
  numerical convergence.

## Current revision status update — 2026-06-19

The build and visual baseline above describes the imported draft at the start
of the revision pass.  The active revision is
`Paper/revision_2026-06-19/draft.tex`; its current status is recorded in
`progress_report.md` and `integrity_audit_2026-06-19.md`.

Key updates since the initial inventory:

- The manuscript now includes an abstract, conclusion/limitations, data and
  code availability, declarations, AI-assisted-preparation disclosure, and a
  numerical reproducibility summary.
- The long-chain placeholder figures were replaced by source-traced figures:
  `action_profiles.pdf`, `cascade_embedding.pdf`, and
  `lte_residual_midchain.pdf`.
- The corrected current-scaling workflow is documented under
  `experiments/flux_validation/`, with canonical simulator/source hashes,
  validation gates, production data, and the manuscript scaling figure.
- LTE table values, short-chain diagnostics, eigen-relaxation diagnostics, and
  manuscript figure metrics are recorded in machine-readable JSON artifacts and
  supporting scripts.
- References were converted to `references.bib` and checked for dangling or
  orphan citation keys.
- The local core numerical claim audit passes 21/21 registered checks; see
  `manuscript_claim_audit.md`.
- The 2026-06-20 submission-level pass further added the Monte Carlo validation
  protocol, the LTE residual mesh figure, timestep sensitivity table, and
  finite-window current diagnostics.  The current one-command gate reports
  `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION` with
  38/38 availability paths present, 21/21 registered numerical claims verified,
  9 author/external submission items pending, and the larger-chain,
  bath-temperature, and thermostat-coupling robustness updates included.

The remaining items before formal journal submission are author/journal
confirmations rather than missing local data files: final declarations,
professional originality screening, target-journal formatting, and a final
post-edit compile/audit pass.
