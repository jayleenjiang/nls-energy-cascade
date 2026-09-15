# R11 fine three-moment calibration report

## Verdict

PASS.  The fixed three-moment exponential calibration passes fine-timestep
validation, the declared final fine test, the equilibrium known-answer
control, and the three-moment fine/coarse timestep gate.  No neural weight,
trajectory, split, shrink factor, bin, diagnostic, or threshold was changed.

## Test access

`RECOVERY_PROTOCOL_R11_FINE_THREE_MOMENT.md` declared before validation that a
second and final V6 test read would be allowed only after both R11 validation
controls passed.  The driven and equilibrium validation verdicts passed.  The
frozen validation hashes were then recorded in
`FORMAL_BLIND_TEST_R11_FINE.md`, after which
`run_recovery_r11_fine_test.sh` read the test once.  No retry occurred.

## Driven test diagnostics

All raw rows are in
`recovery_r11_fine/analysis/diagnostic_comparison.csv`.

| protocol | base seed | ESS | FP median | FP p90 | max moment abs-z | max TV |
|---|---:|---:|---:|---:|---:|---:|
| fine, two moments | 8201 | 0.916008 | 0.088747 | 0.406029 | 1.342781 | 0.024294 |
| fine, two moments | 8202 | 0.916508 | 0.095015 | 0.408698 | 1.443345 | 0.024401 |
| fine, two moments | 8203 | 0.917796 | 0.091415 | 0.396536 | 1.220313 | 0.024247 |
| fine, three moments | 8201 | 0.915647 | 0.088694 | 0.400716 | 1.368753 | 0.024319 |
| fine, three moments | 8202 | 0.916244 | 0.095482 | 0.413261 | 1.450155 | 0.024421 |
| fine, three moments | 8203 | 0.917330 | 0.091638 | 0.401961 | 1.232539 | 0.024276 |
| coarse, three moments | 8201 | 0.908708 | 0.091214 | 0.408893 | 2.829839 | 0.024090 |
| coarse, three moments | 8202 | 0.909603 | 0.089429 | 0.402181 | 2.694922 | 0.024052 |
| coarse, three moments | 8203 | 0.911359 | 0.093102 | 0.434454 | 2.651876 | 0.024189 |

The largest fine three-moment pairwise centered log-ratio RMS is 0.012645,
compared with 0.012733 for the fine two-moment result and 0.013834 for coarse
three-moment replication.  Every per-model and ensemble gate passes.

## Timestep comparison under one protocol

The three-moment fine/coarse ensemble centered log-density RMS is 0.024355 on
the fine test support and 0.024353 on the coarse test support.  Both pass the
fixed 0.10 gate.  The previous asymmetric two-versus-three-moment values were
0.024257 and 0.024272.  Thus applying the third moment at the fine timestep
moves the statistics by only 0.000098 and 0.000082, respectively.

## Calibration coefficients and density shift

The complete standardized and physical coefficients are in
`recovery_r11_fine/analysis/calibration_coefficients.csv`.  The shrunk
standardized coefficients are:

| base seed | two-moment alpha | three-moment alpha |
|---:|---|---|
| 8201 | (-0.012283, 0.002684) | (-0.012219, 0.002679, -0.001261) |
| 8202 | (-0.007364, 0.003292) | (-0.007317, 0.003289, -0.000940) |
| 8203 | (-0.010441, 0.005026) | (-0.010357, 0.005023, -0.001690) |

The first two coefficients barely move.  On 50,000 deterministic fine test
states, the three-minus-two-moment ensemble centered log-density shift has RMS
0.001518, absolute-shift p90 0.002103, p99 0.006155, and maximum 0.023228.
After applying the separately estimated ratio normalizers, the RMS log-density
shift is 0.001518 and its mean is -2.16e-6.  The corresponding centered RMS on
the coarse test support is 0.001543.  Raw per-seed values are in
`recovery_r11_fine/analysis/density_shift.csv`.

## Equilibrium known answer

The exact zero correction was represented by the same three-moment model
class with zero base weights and zero calibration coefficients.  On the final
test, ESS is 1, FP median is 8.88e-16, FP p90 is 3.55e-15, maximum moment
abs-z is 1.272210, maximum TV is 0.033212, and seed disagreement is exactly
zero.  This validates the operator implementation, diagnostic pipeline, and
sampled Gibbs control.  Because the correction is fixed to zero by
construction, it does not test the learned nonequilibrium correction.

## Claim update

The earlier fine/coarse density claim is unchanged in magnitude but is now
procedurally stronger: both timesteps use the identical three-moment,
0.625-shrunk calibration protocol.  The numerical NESS remains validated on
sampled stationary support; no claim is made for unsampled extreme tails or a
global analytic solution.

