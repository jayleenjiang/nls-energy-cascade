# Amendment 003: selector archive filtering

Date: 2026-09-14

The complete repaired 18-fit equilibrium matrix finished successfully.  The
frozen selector then stopped before scientific selection because its glob also
read the two directories preserved by Amendment 002 with the suffix
`_failed_serialization_20260914`.  This produced duplicate records for seeds
4101 and 4102 in one architecture and triggered the exact seed-matrix audit.
The failed selector log is preserved as
`stationary_density/equilibrium_candidates/selection_failed_archive_filter.log`.

The analysis-only correction excludes directories whose name contains
`_failed_serialization_` before grouping metrics.  It does not modify any
candidate metric, seed, architecture, validation observation, admissibility
gate, tie rule, or result.  The required matrix remains the 18 repaired
directories containing both `validation_metrics.json` and `model.weights.h5`.

