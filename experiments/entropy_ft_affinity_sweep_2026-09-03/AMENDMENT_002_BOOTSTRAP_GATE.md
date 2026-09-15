# Analysis amendment 002: bootstrap gate re-specification

Date: 2026-09-04

This amendment is analysis-only.  It was requested after the original
all-resolved-time joint bootstrap returned 0/1000 accepted intercepts for every
case, including the exactly known `(6,6)` equilibrium control.  No trajectory
was regenerated and no production parameter, raw datum, histogram spacing
(`dx = std/20`), minimum count (`10`), minimum number of symmetric pairs (`3`),
or straddle-zero window-selection rule was changed.

## Old gate retained for audit

The old calculation jointly bootstraps every time whose full-sample matched-bin
fit resolves.  Its output and acceptance count remain in the revised tables.
The complete pre-amendment analysis, figures, report, validation note, and
manifest are preserved under `pre_gate_respec_2026-09-04/`.

## New per-time bootstrap

For each `(case,t)`, the same 128 stream resamples are histogrammed on that
time's frozen full-sample bin grid.  The unchanged contiguous reliable-window,
straddle-zero, reflection, `minCount=10`, and `minPairs=3` rules are applied to
that replicate alone.  A replicate is accepted exactly when that time yields a
valid matched-bin slope.  The 2.5 and 97.5 percentiles of accepted slopes form
the reported per-time interval.  The accepted count and whether the interval
contains the case's `Delta beta` are reported for every full-sample-resolved
time.

## New reliable-time intercept bootstrap

For non-equilibrium cases, the time set is frozen from full-sample diagnostics:

```
n_negative >= 500 and matched-bin R^2 >= 0.98.
```

The full-sample slope is regressed on `1/t` over that set.  A stream-bootstrap
replicate is accepted when every time in that set independently resolves under
the unchanged matched-bin rule.  At least 800/1000 accepted replicates and at
least three reliable times are required for a formal intercept verdict.

At equilibrium, the literal `R^2 >= 0.98` set is also reported, but it is empty:
the exact reference relation is a flat zero response, for which regression
`R^2` is not expected to approach one.  To make the known-answer test
well-posed, the equilibrium intercept uses all full-sample-resolved times with
`n_negative >= 500`, omitting only the inapplicable nonzero-signal `R^2` gate.
This exception is explicit in the output and does not apply to driven cases.

## Claim rule

For a driven case, the revised result is consistent with the heat-flux FT only
if at least three reliable times exist, at least 800/1000 joint replicates are
accepted, and the percentile interval for `a_inf / Delta beta` contains one.
At equilibrium the analogous interval must contain zero.  Otherwise the result
is `FAIL` or `UNRESOLVED` according to the same conditions; no point estimate is
promoted without its gate.
