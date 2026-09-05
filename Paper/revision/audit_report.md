# Initial correctness and rigor audit

Date: 2026-06-19

Status: blocking issues found. The current draft is a useful research draft,
but it is not yet suitable for journal submission.

## P0: model/data consistency blockers

1. **The June 18 flux data use the wrong stochastic model.**
   `flux/NLS_flux_SIMD.cpp` omits the required `sqrt(2)` factor in both
   boundary noise amplitudes. This is the same previously documented bug
   fixed in `cpp/simulation/main_fixed.cpp` and in the canonical LTE/5D code.
   Consequently, the reported `n=10,20,30,40` current scaling and tail fits
   cannot be cited as results for the Gibbs-preserving SDE in the manuscript.

2. **The June 18 confidence intervals are computed incorrectly.**
   The code uses `1.96 * std_dev` instead of
   `1.96 * std_dev / sqrt(n_traj)`. The CSV intervals therefore describe
   neither a confidence interval for the mean nor the sample distribution in
   a clearly labeled way.

3. **The boundary phase drift in the manuscript has the wrong sign for its
   stated angle convention.**
   The natural-bath derivation in `Summary_NLS.pdf` and the corrected code use
   a drift proportional to
   `-2 gamma I_neighbor sin(theta)` when
   `theta = 2(phi_neighbor - phi_boundary)`.

4. **The Fokker--Planck diffusion discussion mixes full phase coordinates and
   reduced bond-angle coordinates.**
   The coefficient `4 gamma T/I` is appropriate for a doubled boundary bond
   angle, not for an individual `phi_j` as the surrounding text currently
   states. The state space and diffusion matrix must be written separately in
   the two coordinate systems.

## P0: contradictory numerical claims

1. **Spectral gap:** the manuscript gives approximately `-0.93`, its figure
   caption gives approximately `-0.17`, and the April report gives
   approximately `-1.4`. These values come from different fits and cannot all
   represent the same spectral gap.

2. **“Real spectrum” / overdamping:** a nonoscillatory fit for one or two
   observables does not prove that the generator spectrum is real. At most,
   the measured observable response is compatible with a dominant real decay
   over the fitted window.

3. **Symmetry breaking:** the text reports maximum/mean asymmetry of about
   `67%/28%`, while the figure caption reports about `15%`. The metric,
   mask/normalization, action slice, and value must be recomputed from the
   source array.

4. **Flux scaling:** the draft still states `J ~ n^-1.75` from data generated
   with an older sign-bugged integrator. The June 18 replacement states
   `J ~ n^-1.91` but was generated with the missing-noise-factor bug. Neither
   exponent is currently publication-ready.

## P1: mathematical/statistical issues

1. The manuscript repeatedly calls `J_j` an energy or heat current, while the
   continuity equation shown is for the action/mass `I_j`. The terminology
   must be corrected or an actual Hamiltonian-energy current must be derived.

2. The Fourier-law sign convention is inconsistent: with `T_1 > T_n` and the
   stated positive current direction, defining the gradient as
   `(T_1-T_n)/n` conflicts with a limiting ratio equal to `-kappa`.

3. The Gallavotti--Cohen formula is written using probabilities of exact values
   of a continuous random variable and is asserted without checking the
   hypotheses for this degenerate diffusion. This subsection currently has no
   data and should not make a theorem-level claim about the simulated model.

4. The phase-locking formula uses the wrong inverse-sine branch. The displayed
   expression tends to `pi/3` as `gamma -> 0`, whereas the stable point tends
   to `2pi/3`. The stable branch is
   `pi - asin((sqrt(4 gamma^2+3)-gamma)/(2(1+gamma^2)))`.

5. The flux “tail” data are distributions of per-trajectory finite-time
   averages over a window of length 200, not instantaneous currents. Any
   large-deviation statement must include the averaging window and test
   multiple window lengths.

6. Histogram-bin counts are strongly time-correlated. Poisson/multinomial
   errors based on the raw count total substantially understate uncertainty
   unless effective sample sizes or trajectory/block bootstraps are used.

7. Claims that LTE improves with chain length rely mainly on increasing
   weighted `R^2`. This is not enough to establish convergence to an
   equilibrium marginal; a distance metric with uncertainty and a fixed
   comparison domain is needed.

8. Fast trigonometric/logarithm approximations, positivity clipping, and
   timestep dependence need a systematic bias/convergence study. The existing
   `n=100` artifacts already show that discretization error is material.

## P1: writing and completeness

- No abstract or conclusion.
- No consolidated numerical-methods table.
- No random-seed/reproducibility protocol.
- No data/code availability statement.
- No limitations section.
- Long-chain profile, frequency embedding, LTE residual, and current-scaling
  figures are placeholders.
- Entropy production/fluctuation material is proposed but not implemented.
- The paper does not clearly separate the natural Gibbs-preserving thermostat
  studied numerically from the proof-oriented thermostat in the HLNS paper.
- Several “for the first time” and convergence claims are stronger than the
  current evidence supports.

## Citation metadata verified on 2026-06-19

- HLNS remains arXiv:2505.16018 **v1**, submitted 2025-05-21.
- Zhai--Dobson--Li is PMLR 145, pages 568--597, **2022**.
- CKSTT is Inventiones Mathematicae 181, pages 39--113 (2010),
  DOI 10.1007/s00222-010-0242-2.
- Gallavotti--Cohen is Physical Review Letters 74, 2694--2697 (1995),
  DOI 10.1103/PhysRevLett.74.2694.

## Required evidence before submission

1. Regenerate flux scaling with the corrected canonical SDE, correct standard
   errors, fixed reproducible seeds, and at least one timestep check.
2. Reanalyze the finite-time current distribution over multiple averaging
   windows; avoid calling it an instantaneous-current tail.
3. Resolve the spectral-gap estimator by held-out prediction and fit-window
   sensitivity, or demote the section to exploratory results.
4. Recompute every number appearing in the short-chain figures from saved
   arrays and record the analysis scripts.
5. Add trajectory/block-bootstrap uncertainty for LTE and transport metrics.
6. Replace all figure placeholders and add numerical-convergence appendices.

## Resolution update — 2026-06-19

This file is the initial audit and is retained as a historical record.  The
current manuscript state is tracked in `progress_report.md`,
`integrity_audit_2026-06-19.md`, and
`submission_readiness_checklist_2026-06-19.md`.

Initial blocking issues have been addressed as follows:

| Initial issue | Current resolution status |
|---|---|
| Invalid June 18 flux workflow | Replaced by the canonical Gibbs-preserving simulator and production validation artifacts. The manuscript now reports the corrected finite-size action-current scaling. |
| Incorrect confidence intervals | Replaced by trajectory-standard-error and bootstrap exponent summaries in the canonical workflow. |
| Boundary phase-drift sign and diffusion convention | Corrected in the manuscript's SDE/specification narrative and claim audit. |
| Heat-current terminology | Standardized to action/mass current unless a Hamiltonian-energy current is explicitly discussed as future work. |
| Unsupported entropy/fluctuation-theorem claims | Removed from the results and retained only as a limitation/future-work item. |
| Phase-locking branch | Replaced by the stable reduced fixed-point branch. |
| Spectral-gap/eigenfunction overclaiming | Demoted to an observable-dependent relaxation diagnostic with window sensitivity. |
| Short-chain quantitative inconsistencies | Unsupported percentages and angular-width claims were removed; saved-model diagnostics were rerun and source-traced. |
| Missing abstract, conclusion, availability, limitations, and reproducibility material | Added to `draft.tex`, with a compact numerical reproducibility summary and supporting audit artifacts. |
| Figure placeholders | Replaced by source-traced figures or removed from the manuscript claim set. |

Remaining pre-submission blockers are external or author-supplied:
funding/competing-interest/contribution confirmations, target-journal
formatting decisions, professional plagiarism/self-plagiarism screening, and a
final audit after those edits.
