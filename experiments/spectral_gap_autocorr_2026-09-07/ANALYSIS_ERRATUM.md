# Analysis-only erratum

The first analysis completed successfully and its outputs are preserved under
`analysis_initial_v1/`.  Inspection identified one diagnostic issue: because
the damped-fit frequency is constrained to be nonnegative, a percentile
bootstrap interval can have a tiny positive lower endpoint even when the
optimum is numerically indistinguishable from zero.  Calling that a resolved
oscillation is invalid.

The correction does not alter simulation data, lag grids, autocorrelations,
plateau selection, or plateau rates.  A nonzero frequency is now labelled
resolved only when all three non-tuned identifiability conditions hold:

1. the bootstrap lower endpoint is positive;
2. the fitted interval contains at least half a cycle,
   `lambda_I >= pi/(t_end-t_start)`;
3. the damped model improves AIC by at least 10 relative to its nested
   zero-frequency exponential model.

The rerun also adds the requested descriptive last-supported local rate and
its five-point terminal drift.  These columns are not used to select or alter
the frozen plateau.
