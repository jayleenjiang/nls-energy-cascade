# R11 frozen protocol: identical three-moment calibration at the fine timestep

Frozen before any R11 validation or test statistic is computed.

## Purpose

R11 removes the procedural asymmetry in the published timestep comparison.
The existing fine neural weights are not retrained.  Each fine base model is
wrapped in the same three-observable exponential calibration used by the R10
coarse recovery:

1. `I1`;
2. `I1*I2*sin(theta1)`;
3. `I2*I3*sin(theta3)`.

The fitted calibration vector is multiplied by the already inherited common
shrink factor `0.625`.  No alternative factor, observable, model seed, split,
normalizer, diagnostic, or gate will be examined.

## Frozen inputs

- base neural seeds: `8201, 8202, 8203`;
- reported R11 seeds: `47201, 47202, 47203`;
- base neural weight SHA-256 values:
  - `7c76f5cbd00b01db27a1d9f831c4b809c4db95deff47079a43a69a9f937ea322`;
  - `881b96ef4aa69cc94ac6f8b7ca239cbbefdcc8a4d582f0efa0802587aa221738`;
  - `8fd45218b837b7026b6128d7e2578f06415ce1ca4986b0bf6faf392dc447084f`;
- fine holdout manifest: `V6_HOLDOUT_SPLITS.csv`, SHA-256
  `b08bb058562b108cf674cad2dff2366bfafa9d7f98b1ca7d03f009bef668f914`;
- coarse holdout manifest: `V4_HOLDOUT_SPLITS.csv`, SHA-256
  `1c51aa91d15f076e7b411152425fd6dbd6ab34d28f9b473d7da10409d80cc2c5`.

The calibration estimation rows are exactly those used by R8/R9: all already
opened legacy fine train/validation/test streams and V4 fine
validation/test streams, separately for the driven target and equilibrium
reference.  The V6 validation and test streams are not used to fit the
calibration.

## Unchanged gates

Each of the three models must independently satisfy:

- finite, unclipped importance weights;
- ESS fraction at least `0.10`;
- median and 90th-percentile absolute relative stationary residual at most
  `0.10` and `0.50`;
- all twelve whole-stream-bootstrap moment discrepancies satisfy
  `|z| <= 3.748286557303567`;
- all one- and two-angle marginal total-variation distances are at most
  `0.05`;
- largest pairwise centered log-ratio RMS across seeds is at most `0.10`.

Marginal edges are created on R11 validation and frozen for R11 test.

## Test-access declaration

The V6 fine test split was previously opened once for the two-moment R9 model.
R11 authorizes one additional and final access for the fully specified
three-moment transformation.  This access is conditional on complete R11
driven and equilibrium validation PASS.  Test values cannot alter the model,
calibration, shrink factor, bins, gates, or reporting.  A failure is final for
this protocol and will be retained.

## Equilibrium known answer

The equilibrium control is passed through the same three-moment model class,
with the base neural correction and all three calibration coefficients set
exactly to zero.  It tests the operator and diagnostics and whether sampled
equilibrium marginals match the exact Gibbs reference.  It does not test the
learned nonequilibrium correction.

## Interpretation

The fine/coarse comparison is upgraded only if R11 validation and the declared
test access pass.  The primary timestep statistic is the ensemble centered
log-density RMS on both fine and coarse test supports, with the existing
threshold `0.10`.  The previous two-moment fine artifacts are retained.

