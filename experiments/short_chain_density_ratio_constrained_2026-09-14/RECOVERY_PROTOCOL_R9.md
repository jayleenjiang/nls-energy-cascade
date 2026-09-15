# Recovery protocol R9: predeclared calibration shrinkage

Frozen before R9 candidate output exists on 2026-09-14.

The V6 driven blind test passed every density and moment check except one:
seed 43202 had relative stationary-FP median 0.100997 versus the frozen 0.10
limit.  The other two seeds had 0.096266 and 0.099129.  This is consistent with
the R8 two-moment exponential correction being slightly too large for the FP
margin.  All V5 rows are now opened development data.

The equilibrium exact-zero ratio was rejected by the I3 moment in the V6 test.
An empirical null calibration using 100,000 predeclared 16+16 splits of all
192 opened equilibrium fine streams found `P(max |z| >= 3.7482865573)=0.00707`.
Thus the old cutoff is already conservative; it is retained unchanged.  The
V6 equilibrium result is recorded as a rare finite-sample false rejection and
is not used to loosen the gate.

## Candidate family and selection

Starting from each exact R8 bundle, multiply both frozen calibration
coefficients by one common factor in the fixed grid

    f in {0.500, 0.625, 0.750, 0.875, 1.000}.

No neural weight, center, scale, feature, or relative coefficient is changed.
Evaluate every candidate on the four opened development splits: V4 validation,
V4 test, V5 validation, and V5 test.  Use the unchanged formal diagnostics and
gates.

For method selection only, define the safety score as the maximum over all
splits and models of

    FP_median/0.09, FP_p90/0.45, max_abs_moment_z/3.748286557303567,
    max_marginal_TV/0.05, pairwise_log_ratio_RMS/0.10.

Select the factor with the smallest safety score among candidates that pass
all original formal gates on all four opened splits.  Ties are resolved in
favor of the smaller factor.  The safety denominators do not replace the
formal gates and will not be changed after seeing results.

If no candidate passes all original gates, R9 fails.  If R9 succeeds, generate
a fresh V6 fine heldout dataset before any new formal verdict.
