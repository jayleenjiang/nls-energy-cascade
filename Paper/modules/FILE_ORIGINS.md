# Canonical file origins and duplicate-directory policy

The repository contains historical workspaces, curated mirrors, production
experiment directories, and manuscript-ready outputs.  Similar filenames do
not necessarily mean interchangeable evidence.

| Path | Role | Use policy |
|---|---|---|
| cpp/ | Canonical C++ source grouped by simulation, fp5d, and backward tasks | Prefer for current source citations |
| python/ | Canonical Python analysis, plotting, and data-generation source | Prefer for reusable analysis |
| NN notebooks/ | Readable notebook copies for neural Fokker--Planck and eigenfunction work | Prefer when explaining algorithms |
| KDE/ | Historical short-chain working directory containing raw outputs, saved models, and notebook results | Primary source for archived short-chain numerical artifacts |
| experiments/kde/ | Curated mirror of many KDE experiments | Use for browsing; check hashes before replacing KDE sources |
| lte/ | Primary long-chain LTE histograms, marginals, profiles, and conditioned data | Primary source for LTE tables and residuals |
| experiments/lte/ | Curated profile/analysis mirror with lighter outputs | Use for scripts and manuscript profile generation |
| flux/ | Current action-current simulators and analysis scripts | Prefer NLS_flux_canonical.cpp for scaling and NLS_flux_relaxation_tau.cpp for burn-in/tail work |
| flux_data/ | Historical canonical/SIMD comparison datasets | Comparison only; not the primary production location |
| Paper/revision/experiments/flux_validation/ | Dated, parameterized production and validation results | Primary source for current current-scaling and burn-in studies |
| Paper/revision/ | Active locally enhanced manuscript, figures, metrics, and audits | Current local paper workspace |
| Paper/packages/ | Frozen compile/share packages | Do not edit in place |
| Paper/modules/ | This section-oriented navigation and provenance layer | Update whenever claims or experiment paths change |

## Current versus historical revision paths

Older reports and JSON files often contain Paper/revision_2026-06-19.  Git
still records deletions at that historical path, while the live local files
are under Paper/revision.  Do not restore, delete, or bulk-stage either tree
without a separate migration decision.

## Flux code policy

- Primary scaling: flux/NLS_flux_canonical.cpp.
- Latest burn-in and finite-time distribution: flux/NLS_flux_relaxation_tau.cpp.
- SIMD comparison: flux/flux_data/NLS_flux_SIMD.cpp and dated fixed-SIMD data.
- Older cpp/simulation/flux_V1.cpp and flux_V2.cpp are legacy references and
  must not be used to support current paper numbers.

## LTE MATLAB source

The advisor-requested residual mesh was produced from the MATLAB workspace at
/Users/jayleenjiang/Documents/MATLAB/lte.  Preserved copies of
compare_residual.m, residual.m, and fit_new.m are stored in
02_long_chain_ness_lte/code/matlab so the paper module is self-contained.
