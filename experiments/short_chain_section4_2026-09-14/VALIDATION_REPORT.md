# Section 4 validation report

## Frozen design

The analysis was registered in `PROTOCOL.md` before model fitting. Physical
trajectories were partitioned by trajectory identity, never by shuffled
snapshots: 32 training, 16 validation, and 16 blind-test streams per bath and
timestep. Density and mode branches used separate deterministic split salts.
The primary fine-timestep blind-test effective sample sizes were 11,616
(driven) and 13,860 (equilibrium), above the frozen minimum of 5,000.

The reduced generator audit passed: maximum drift error
`2.27e-13`, maximum covariance error `2.91e-11`, equilibrium Gibbs relative
adjoint residual RMS `1.51e-14`, and periodic-energy error `8.53e-14`.
The direct-state replay preserved all 256 source trajectories with maximum
observable discrepancy `1.14e-13` and removed a float32 reconstruction
cancellation without altering the dynamics.

## Equilibrium weak-identity calibration

The initial test family contained six bare angular functions for which the
generator terms are not square-integrable near zero action. This was detected
before density fitting and is preserved in
`ANALYSIS_ERRATUM_001_WEAK_DOMAIN.md`. The frozen simultaneous calibration was
then applied to the remaining 18 regular functions. Its family-wise threshold
was 3.748 standard errors. No equilibrium or driven blind direct-trajectory
identity failed. This validates the trajectory/operator consistency check; it
does not rescue a density that fails the pointwise adjoint equation.

## Equilibrium density candidate matrix

All entries below summarize three independent seeds. The selection gates were
Gibbs slope in `[0.95,1.05]`, median pointwise residual at most `0.20`, and
90th-percentile residual at most `1.00` for every seed.

| K | widths | lambda_FP | median val. NLL | Gibbs slope range | centered log-RMSE range | median residual range | p90 residual range | admissible |
|---:|:---:|---:|---:|:---:|:---:|:---:|:---:|:---:|
| 8 | 64x64 | 0.01 | 7.9010 | [1.0129,1.0370] | [0.4558,0.5471] | [1.886,2.624] | [8.096,9.456] | no |
| 8 | 64x64 | 0.10 | 7.9487 | [1.0117,1.0474] | [0.5454,0.5954] | [2.018,2.478] | [8.044,9.277] | no |
| 16 | 64x64 | 0.01 | 7.8255 | [1.0268,1.0592] | [0.3034,0.3520] | [1.253,1.419] | [5.075,6.075] | no |
| 16 | 64x64 | 0.10 | 7.9112 | [1.0126,1.0401] | [0.4646,0.5319] | [1.133,1.610] | [3.705,6.085] | no |
| 16 | 128x128 | 0.01 | 7.8291 | [1.0342,1.0696] | [0.3020,0.3533] | [1.174,1.406] | [4.552,5.698] | no |
| 16 | 128x128 | 0.10 | 7.9494 | [1.0162,1.0628] | [0.4511,0.5860] | [1.285,1.520] | [3.772,5.529] | no |

The pointwise residual misses the frozen gates by factors of at least 5.7
(median) and 3.7 (p90). The outcome is therefore robust to choosing the most
favorable architecture or seed. Candidate selection stopped with
`NO_ADMISSIBLE_CANDIDATE`; blind equilibrium density evaluation and NESS
density training were not performed.

## Held-out generator-mode analysis

The mode branch proceeded independently of the density failure. The slow real
mode used the frozen long-lag selection; the oscillatory mode used the frozen
short-lag selection and frequency band. Whole-trajectory bootstrap used 500
replicates.

| case | selected EDMD Re(lambda) | blind Re(lambda) [95% CI] | blind Im(lambda) [95% CI] | test residual | validation residual |
|:---|---:|:---:|:---:|---:|---:|
| driven, dt=1e-3, real | -0.9139 | -0.8840 [-0.9707,-0.8103] | -- | 0.769 | 0.781 |
| driven, dt=2.5e-4, real | -0.9515 | -0.9088 [-1.0139,-0.8178] | -- | 0.784 | 0.781 |
| equilibrium, dt=1e-3, real | -0.9376 | -0.9507 [-1.0751,-0.8420] | -- | 0.774 | 0.776 |
| equilibrium, dt=2.5e-4, real | -0.9086 | -0.8264 [-0.9083,-0.7587] | -- | 0.773 | 0.777 |
| driven, dt=1e-3, complex | -4.7880 | -7.067 [-10.804,-4.452] | 5.059 [0.0003,6.983] | 0.807 | 0.850 |
| driven, dt=2.5e-4, complex | -5.1704 | -7.198 [-8.885,-6.401] | 6.281 [5.388,6.862] | 0.668 | 0.733 |
| equilibrium, dt=1e-3, complex | -5.6755 | -6.986 [-11.137,-5.563] | 4.438 [0.024,6.604] | 0.751 | 0.794 |
| equilibrium, dt=2.5e-4, complex | -5.4079 | -5.967 [-7.085,-5.261] | 5.088 [4.487,5.558] | 0.695 | 0.713 |

The blind complex fits strongly prefer the damped model by AIC, but the coarse
frequency intervals are too broad for precise numerical comparison. The
fine-timestep fits agree with the independent observable-level frequencies.

## Audit trail

- `AMENDMENT_001_DIRECT_STATE_OUTPUT.md`: float64 direct-state replay after a
  saved-representation cancellation, before model fitting.
- `AMENDMENT_002_MODEL_SERIALIZATION.md`: serialization-only repair after two
  completed fits failed while saving; trainable arrays were unchanged.
- `AMENDMENT_003_SELECTOR_ARCHIVE_FILTER.md`: selector ignores preserved
  failed-serialization directories; no scientific metric changed.
- `ANALYSIS_ERRATUM_001_WEAK_DOMAIN.md`: excludes six infinite-variance weak
  probes before density fitting.
- `ANALYSIS_ERRATUM_002_MODAL_WEIGHTS.md`: deterministic held-out modal-weight
  read added without reselection or gate changes.

The complete per-seed density table is
`analysis/density_candidate_metrics.csv`; held-out mode values are in
`analysis/heldout_mode_summary.csv`; modal weights are in
`modes/blind_modal_weights.csv`.

