# Revision roadmap

## Gate A: canonical model and reproducible simulation

- Write one authoritative SDE specification in both full and reduced
  coordinates.
- Add a unit/sanity test showing equilibrium agreement with
  `exp(-H/(2T))`.
- Repair the flux simulator and record deterministic seeds and run metadata.
- Run short equilibrium and NESS smoke tests before expensive production runs.

## Gate B: long-chain evidence

- Produce converged action profiles for multiple chain lengths.
- Recompute current scaling with valid uncertainty and timestep sensitivity.
- Analyze finite-time averaged-current distributions at several window lengths.
- Recompute LTE distances with block/bootstrap uncertainty and fixed masks.
- Decide whether entropy production is completed or removed from the paper.

## Gate C: short-chain evidence

- Audit the 5D neural Fokker--Planck loss, normalization, validation slices,
  and held-out error metrics.
- Recompute symmetry-breaking and phase-locking statistics.
- Resolve or substantially narrow the generator/eigenfunction claims.

## Gate D: manuscript

- Add abstract, conclusion, limitations, methods table, reproducibility
  statement, and data/code availability.
- Replace every placeholder with a source-traceable figure.
- Move exploratory or weakly supported claims to a clearly labeled discussion.
- Standardize terminology: action/mass current versus Hamiltonian energy.
- Convert references to a maintainable BibTeX database with verified metadata.

## Gate E: final quality

- Compile cleanly with no unresolved references or overfull boxes.
- Visually inspect every PDF page.
- Run a fresh citation/data/claim integrity audit.
- Perform an independent reviewer-style methodology and presentation review.
- Commit and push only coherent, reproducible milestones.

## Gate status update — 2026-06-19

This roadmap is retained as the initial repair plan.  Current status:

| Gate | Status | Notes |
|---|---|---|
| Gate A: canonical model and reproducible simulation | Complete for manuscript claims | Canonical Gibbs-preserving current simulator, validation gates, reproducible seeds/metadata, and source hashes are recorded. |
| Gate B: long-chain evidence | Complete within narrowed claim scope | Action profiles, corrected current scaling, LTE diagnostics, and reproducibility artifacts are included. Entropy-production/GC and open-chain terminal-energy claims are intentionally not claimed. Larger-chain and smaller-timestep studies remain future work. |
| Gate C: short-chain evidence | Complete within saved-model diagnostic scope | Equilibrium validation, qualitative symmetry breaking, phase-locking branch, and eigen-relaxation diagnostics are source-traced. Full neural-network retraining remains optional if requested by the journal. |
| Gate D: manuscript | Complete except author-supplied declarations | Abstract, conclusion, limitations, availability, reproducibility summary, source-traced figures, and BibTeX references are in place. Funding/contribution/competing-interest confirmations still require author review. |
| Gate E: final quality | Partially complete; final pass still required | Current local claim audit passes, and major updates are pushed to GitHub. After author/journal edits, rerun compile, visual PDF QA, citation/data/claim audit, and professional originality screening. |

See `submission_readiness_checklist_2026-06-19.md` for the remaining
author-facing and journal-facing checklist.

## Submission-level update — 2026-06-20

The 2026-06-20 pass strengthened the journal-facing evidence without expanding
the claim scope: the manuscript now includes the Monte Carlo validation and
uncertainty protocol, the `n=15,25,50` LTE residual mesh diagnostic from the
`compare_residual.m` convention, a current-estimator timestep sensitivity
table, and finite-window current diagnostics.  The current local gate reports
`PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`:
LaTeX/log checks pass, availability checks cover 40/40 paths with zero
untracked required files, reference integrity passes, and the registered
numerical claim audit verifies 23/23 claims. The author/submission-field audit
records 9 pending author or external-service items. Remaining blockers are
author/journal/external-release decisions:
target journal/template, final declarations, professional similarity screening,
and whether to create a DOI-backed raw-data archive.

## Larger-chain robustness update — 2026-06-20

An additional `n=50` current run has been completed with the canonical
Gibbs-preserving current accumulator (`1024` trajectories, burn-in `8000`,
measurement window `200`, `dt=5e-4`).  The result supports the insulating
transport conclusion and steepens the diagnostic five-size fit to exponent
`-1.89449`.  A later `dt=2.5e-4` production-resolution check at the same length
gives a compatible mean current (`3.60%`, or `1.12` pooled standard errors,
above the `dt=5e-4` run).  A matched production-resolution `n=60`,
`dt=2.5e-4` check was then completed, giving a compatible `3.50%` upward shift
relative to the `dt=5e-4` `n=60` run (`0.78` pooled standard errors).  The
larger-length results are retained as robustness checks rather than replacing
the primary `n=10,20,30,40` exponent.  The immediate larger-length/timestep
sanity checks are now addressed; a full asymptotic convergence study would
still require still larger lengths and/or a dedicated finite-timestep
convergence campaign.
