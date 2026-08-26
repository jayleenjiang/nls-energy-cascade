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

## Analysis outputs

Source: `../../flux/analyze_entropy_ft.py`.

The `analysis_v2/` directory contains:

- symmetric-bin diagnostics for medium entropy, bath heat, and action current;
- both-tail survival probabilities on an `A=0.01,0.02,...` grid;
- raw and plus-four probabilities;
- fitted-normal survival benchmarks for both action-current tails;
- stream-bootstrap confidence intervals;
- action/heat correlation and residual-variance diagnostics;
- stationarity and sample-size tables.

The pilot source hashes and exact parameters are recorded in
`pilot/pilot_manifest.txt`.  Raw production blocks remain local because their
expected size is several hundred megabytes; source, manifests, summaries, and
curated figures are the GitHub-facing reproducibility artifacts.
