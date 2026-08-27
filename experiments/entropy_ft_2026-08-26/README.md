# Joint entropy/action fluctuation study (2026-08-26)

## Status

The completed `pilot/` run is **provisional**.  It is a validated sample-size
and observable audit, not a final fluctuation-theorem (FT) production result.
No pilot file is overwritten by the production workflow.

## Questions and observables

The experiment keeps two questions separate while obtaining both observables
from the same Cartesian trajectory.

1. Thermodynamic FT diagnostic: record the heat delivered by each bath,
   `Q_left` and `Q_right`, and form

   `Sigma_medium = -Q_left/T1 - Q_right/Tn`.

   For the rate `a = Sigma_medium/t`, the asymptotic steady-state diagnostic is
   `R_t(a) = log[p_t(a)/p_t(-a)]/t`; a unit slope is the reference only after
   finite-time and system-entropy boundary terms are under control.

2. Action-current fluctuations: record the bulk action current `J_t`.  Its
   symmetry slope is an empirical finite-time rate proxy; it is not identified
   a priori with thermodynamic entropy production.  The joint sample also
   measures its coupling to the bath heat current.

The Hamiltonian convention is `E = H_code/2 = H_paper`, so the inverse bath
temperature multiplying physical heat is `1/T`.  Positive `Q_r` means energy
delivered by bath `r` to the chain.

## Sampler and required audit fields

Source: `../../flux/NLS_entropy_ft.cpp`.

- projection-free Cartesian variables `c_j = x_j + i y_j`;
- no `ACTION_FLOOR` or positivity clipping;
- implicit-midpoint Hamiltonian update;
- independent left/right bath substeps;
- heat accumulated as the physical-energy change across the corresponding
  bath-only substep, not as a whole-step energy difference;
- every `t=20` record contains `Q_left`, `Q_right`, `delta_E`,
  `Sigma_medium`, `stream_id`, `block_id`, action current, and the energy
  balance residual `Q_left + Q_right - delta_E`.

## Numerical validation

The self-test checks the Cartesian gradient, orthogonality of Hamiltonian flow
to the energy gradient, and the boundary Laplacian.  All production parameters
use `T1=10`, `Tn=2`, `gamma=0.1`, and `dt=5e-4`.

Matched `n=40` timestep subsets (1,280 `t=20` blocks each) give:

| dt | mean action current | mean entropy rate | balance RMS rate | midpoint failures |
|---:|---:|---:|---:|---:|
| 5e-4 | 0.03113 | 0.20872 | 2.90e-6 | 0 |
| 2.5e-4 | 0.02892 | 0.20746 | 7.45e-7 | 0 |

The current and entropy means overlap within stream-level Monte Carlo
uncertainty.  Halving the timestep reduces the balance residual by about a
factor of four.  Equal-temperature and bath-swapped controls are stored in
`validation/`; the former has zero mean current within uncertainty and the
latter reverses the action current while retaining positive entropy production.

## Completed pilot

Each chain length has 128 independent streams and 20,096 non-overlapping
`t=20` blocks after burn-in `B=500`.

| n | mean action current | mean heat current | mean entropy rate | negative entropy blocks | balance RMS rate |
|---:|---:|---:|---:|---:|---:|
| 10 | 0.40143 | 3.24088 | 1.29683 | 143 | 9.73e-6 |
| 20 | 0.11968 | 1.40515 | 0.56219 | 2,201 | 5.68e-6 |
| 30 | 0.05486 | 0.80402 | 0.32151 | 4,327 | 4.05e-6 |
| 40 | 0.03019 | 0.51839 | 0.20736 | 5,590 | 3.12e-6 |

The mean action currents agree with the established canonical transport data.
Stationarity checks show no material first-to-last-quarter drift.  The pilot
does **not** establish the FT: the directly fitted medium-entropy symmetry
slopes in the resolvable `t=20` windows are below one, and several larger-time
windows do not contain enough negative events for an unsmoothed fit.  This is
reported as a diagnostic result rather than corrected or renormalized away.

## Production decision

The production target is 1,000,064 non-overlapping `t=20` blocks per chain
length (`n=10,20,30,40`).  The count is set by the rarest `n=10` negative tail:
the pilot probability predicts roughly 7,100 negative medium-entropy and
17,000 negative action-current blocks at `t=20`.  Aggregated larger-`t` windows
will be reported only where raw positive and negative bin counts pass the
predeclared effective-count threshold; plus-four probabilities are for plots,
never for the inferential fit.

The pilot also predeclares the expected resolution limit of the production
sample.  Using 50 effective negative events as the minimum for a two-sided
fit, the largest likely resolvable averaging windows are:

| observable | n=10 | n=20 | n=30 | n=40 |
|---|---:|---:|---:|---:|
| medium entropy rate | 40 | 80 | 120 | 200 |
| action-current rate | 40 | 120 | 200 | 200 |

All requested windows through `t=200` are still generated.  A window beyond
this forecast is reported as unresolved if the production run contains too few
raw negative events; plus-four smoothing is not used to convert that censoring
into a fitted result.  Extending those short-chain, long-time tails by direct
sampling would require exponentially more trajectories and is a separate
rare-event calculation rather than a reason to alter the predeclared fit rule.

## Analysis outputs

Source: `../../flux/analyze_entropy_ft.py`.

The `analysis_v2/` directory contains:

- symmetric-bin diagnostics for medium entropy, bath heat, and action current;
- both-tail survival probabilities on an `A=0.01,0.02,...` grid;
- raw `log P` curves and scaled survival-rate proxies `-log(P)/t` for testing
  finite-time large-deviation collapse;
- raw and plus-four probabilities;
- full-sample normal survival benchmarks and a shared-parameter descriptive
  Gaussian fit to both raw action-current tails;
- stream-bootstrap confidence intervals;
- action/heat correlation and residual-variance diagnostics;
- stationarity and sample-size tables.

The predeclared fixed-width result is retained as the primary symmetry
diagnostic.  `../../flux/analyze_ft_adaptive_bins.py` supplies a separate
fit-range robustness check: it uses symmetric equal-width bins, fixes the
range from the 99th percentile of each sign separately, and limits the bin
count by rare-side effective support.  This avoids the discontinuity that
occurs when the overall first percentile crosses zero.  Differences between
the primary and robustness slopes are reported as curvature/fit-window
sensitivity; the estimate closer to the theoretical reference is never
selected post hoc.  The independent analysis auditor recomputes both results
from the raw blocks.

The pilot source hashes and exact parameters are recorded in
`pilot/pilot_manifest.txt`.  Raw production blocks remain local because their
expected size is several hundred megabytes; source, manifests, summaries, and
curated figures are the GitHub-facing reproducibility artifacts.

## Production acceptance protocol

The production directory is not considered a scientific result merely because
the sampler processes exit successfully.  Acceptance requires all of the
following gates, in this order:

1. `production_manifest.txt` contains `completed_utc` and no `failed_utc`.
2. Each of `n10`, `n20`, `n30`, and `n40` contains exactly 1,000,064 finite
   blocks with the declared stream/block ordering.
3. The recorded entropy and energy-balance columns recompute from
   `Q_left`, `Q_right`, and `delta_energy` to numerical tolerance; all midpoint
   solves succeed; the balance RMS rate remains below the declared audit
   threshold.
4. The source hashes in the launch manifest still match the sampler and core
   analyzer used for the run.
5. Symmetry fits use raw nonzero counts on both sides.  Plus-four estimates are
   display-only, and missing rare tails remain missing rather than being
   extrapolated into evidence.  The primary fixed-range result and the
   adaptive-range robustness result are both retained.
6. The final interpretation checks the trend with averaging time, the number of
   negative events, normal-tail residuals, and heat--action residual variance.
   A failed unit-slope check for medium entropy is reported as “not verified in
   the sampled window,” because the NESS system-entropy endpoint term has not
   been reconstructed.

From the repository root, the post-run commands are:

```bash
flux/finalize_entropy_ft_run.sh
```

This fail-fast wrapper creates `production/final_v1/`, refuses to overwrite an
existing finalization, enforces the production source hashes, runs both audits,
builds the supplementary analysis and report, compiles the PDF, and writes
separate SHA-256 lists for the four raw block files and the curated artifacts.
The expanded commands executed by the wrapper are shown below for transparency:

```bash
/opt/homebrew/bin/python3 flux/audit_entropy_ft_run.py \
  experiments/entropy_ft_2026-08-26/production \
  --expected-blocks 1000064 \
  --require-completed \
  --require-hash-match \
  --output-prefix experiments/entropy_ft_2026-08-26/production/audit

/opt/homebrew/bin/python3 flux/plot_entropy_ft_supplement.py \
  experiments/entropy_ft_2026-08-26/production/n*_blocks.csv \
  --analysis-dir experiments/entropy_ft_2026-08-26/production/analysis \
  --output-dir experiments/entropy_ft_2026-08-26/production/final_v1/supplement \
  --taus 20,40,60,80,100,120,140,160,180,200 \
  --threshold-step 0.01 \
  --minimum-raw-count 5 \
  --tail-probability-max 0.01

/opt/homebrew/bin/python3 flux/analyze_ft_adaptive_bins.py \
  experiments/entropy_ft_2026-08-26/production/n*_blocks.csv \
  --output-dir experiments/entropy_ft_2026-08-26/production/final_v1/adaptive \
  --taus 20,40,60,80,100,120,140,160,180,200 \
  --max-bins 60 \
  --min-effective-count 50 \
  --range-quantile 0.99 \
  --bootstrap 1000

/opt/homebrew/bin/python3 flux/audit_entropy_ft_analysis.py \
  experiments/entropy_ft_2026-08-26/production \
  --analysis-dir experiments/entropy_ft_2026-08-26/production/analysis \
  --supplement-dir experiments/entropy_ft_2026-08-26/production/final_v1/supplement \
  --adaptive-dir experiments/entropy_ft_2026-08-26/production/final_v1/adaptive \
  --output-prefix experiments/entropy_ft_2026-08-26/production/final_v1/analysis_audit

/opt/homebrew/bin/python3 flux/build_entropy_ft_report.py \
  --run-dir experiments/entropy_ft_2026-08-26/production \
  --analysis-dir experiments/entropy_ft_2026-08-26/production/analysis \
  --supplement-dir experiments/entropy_ft_2026-08-26/production/final_v1/supplement \
  --adaptive-dir experiments/entropy_ft_2026-08-26/production/final_v1/adaptive \
  --validation-dir experiments/entropy_ft_2026-08-26/validation \
  --output-dir experiments/entropy_ft_2026-08-26/production/final_v1/report \
  --status-label production
```

The generated report must then be compiled and visually inspected before any
number or conclusion is migrated into the paper.
