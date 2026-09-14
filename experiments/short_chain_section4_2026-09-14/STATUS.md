# Execution status

- [x] Experiment directory isolated from the dirty manuscript worktree.
- [x] Scientific protocol and verdict boundaries frozen before model fitting.
  Protocol SHA-256:
  `0513a546ba80ea0bd8abbbb990c93b74bc416e83732701941429bea391ae1b70`.
- [x] Reduced generator/Fokker--Planck operator independently audited:
  **PASS**.
- [x] Normalized-density automatic-differentiation smoke test: **PASS**.
- [x] Historical input/hash audit and trajectory splits completed.
  Three cases pass. `equilibrium_dt2p5e-4` fails the frozen positivity gate
  because two float32 sum/difference rows reconstruct one action as zero.
- [x] Amendment-001 float64 direct-state replay and same-trajectory audit:
  **PASS**. All 256 replayed streams match the old observables, with maximum
  absolute difference `1.1368683772161603e-13`; no model fit has started.
- [x] Split-level effective sample sizes measured: **PASS**. The conservative
  fine-step test ESS is 11,616 for NESS and 13,860 for equilibrium, above the
  frozen threshold of 5,000.
- [x] Equilibrium weak-identity resolution calibrated and blind direct-
  trajectory control applied: **PASS** on the 18 regular generator-domain
  functions. The corrected simultaneous threshold is 3.748 sigma and there
  are no equilibrium or NESS blind-test failures. The initial invalid 16.805-
  sigma artifact is preserved.
- [x] Section 4.1 equilibrium model matrix completed. All 18 fits finished,
  but none of the six architecture/penalty families passed the frozen
  equilibrium-validation gates for all three seeds. The selector therefore
  wrote `NO_ADMISSIBLE_CANDIDATE`. In accordance with the protocol, the blind
  equilibrium density split was not opened and no NESS density was trained.
- [x] Section 4.3 one-shot held-out mode analysis completed. The frozen slow
  real EDMD-visible modes pass the validation-to-test residual-growth and
  coarse/fine interval-overlap gates. The oscillatory fits also pass the
  predeclared reportability gates, although the coarse-step frequency
  intervals are broad. The unique-spectral-gap claim remains forbidden by the
  pre-existing common-autocorrelation-plateau failure.
- [x] Section 4.2 density-based interpretation stopped by the Section-4.1
  known-answer gate. No quantitative NESS-density or stabilization-mechanism
  claim is made from the failed family.
- [x] Publication figures, validation tables, final verdict, and rewritten
  Section 4 TeX generated.
- [x] Standalone PDF visually verified and the complete reproducibility
  package archived with a tested compressed stream.
- [x] Scoped report/code artifacts committed and pushed to
  `codex/paper-journal-revision` (deliverables commit `0c41a1f`).

The frozen equilibrium density candidate matrix has terminated with the
predeclared negative result `NO_ADMISSIBLE_CANDIDATE`. Blind density-test
statistics were not inspected. Section 4.3 used its test streams exactly once
after freezing `MODE_SELECTION.json`; the immutable selection hash and the
completed read marker are recorded under `modes/`.

The historical failure is documented in `INPUT_INTEGRITY_FAILURE.md`. It is a
saved-representation failure, not an integrator floor: the source trajectory
state is double precision and the run metadata reports a strictly positive
minimum action of `1.1709150955926485e-08` with zero projection/floor events.
