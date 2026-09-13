# Final verdict: finite-time and long-chain fluctuation relations

## Bottom line

The completed data support two distinct numerical statements.

1. At `n=2`, a finite-time NESS **total-entropy** relation passes both a
   parametric endpoint-density audit and an independent discrete path-ratio
   audit.
2. At `n=10,20,30,40`, the long-time entropy-current cloning calculations
   are numerically consistent with the Gallavotti--Cohen SCGF relation over
   the resolved complementary tilt pair at each length.  Every frozen
   support, observation-time, population, timestep, and selection-interval
   gate now passes.

This is a controlled numerical verification in the tested windows.  It is
not a mathematical proof, a test of the full tilt interval, or an
infinite-chain theorem.

## Long-chain result

| n | tested pair | final `N_c,t` | `D=psi(k)-psi(1-k)` | 95% CI | complete controls |
|---:|:---:|:---:|---:|:---:|:---:|
| 10 | (0.3,0.7) | 2048, 60 | -0.000075 | [-0.009953, 0.009802] | pass |
| 20 | (0.4,0.6) | 1024, 60 | -0.003389 | [-0.012081, 0.005302] | pass |
| 30 | (0.4,0.6) | 4096, 60 | -0.002336 | [-0.007251, 0.002578] | pass |
| 40 | (0.4,0.6) | 1024, 120 | 0.004498 | [-0.002343, 0.011340] | pass |

The admissible paper-level statement is:

> For the boundary-driven projection-free Cartesian NLS chain at
> `(T_L,T_R,gamma)=(10,2,0.1)`, controlled rare-event estimates are
> consistent with the Gallavotti--Cohen SCGF symmetry at all four tested
> chain lengths over the resolved complementary tilt pair.  The conclusion
> is stable under the prospectively specified population, observation-time,
> timestep, and selection-interval checks.

The secondary signed residuals normalized by the mean SCGF magnitude are
approximately -0.04%, -4.0%, -4.9%, and 12.9% for `n=10,20,30,40`.  The
`n=40` absolute gate passes, but its relative uncertainty is visibly larger
because the SCGF magnitude is smaller; that caveat must accompany the result.

## Relation to the mentor's direct-sampling question

The requested direct experiment is complete at `n=10,20,30,40`, with
1,000,064 non-overlapping `t=20` blocks per length aggregated to
`t=20,40,...,200`.  It establishes resolved finite-time large-deviation
scaling in the action-current tails.  The directly sampled medium-entropy and
heat symmetry slopes remain below their asymptotic references before the
negative tails disappear.  Therefore direct sampling alone gives the result
"not verified in the accessible finite-time window," not a violation.

The cloning calculation answers the complementary rare-event question: after
all convergence controls, the resolved long-time SCGF pairs are consistent
with Gallavotti--Cohen symmetry.  These statements are compatible because
the direct estimator is endpoint- and rare-tail-limited at long averaging
times.

## Small-chain finite-time control

For `n=2`, the accepted parametric NESS total-entropy analysis used 1,048,576
blocks per condition.  It gives a detailed-FT slope
`0.994496 +/- 0.00738`, a driven integral-FT log estimate `0.002175` with
95% CI `[-0.000744,0.005251]`, and a validated equilibrium density model.
The independent discrete path-ratio audit gives a Crooks slope
`0.993352 +/- 0.00243`; its forward and reverse integral relations both
pass.  This is a finite-time `n=2` verification, not the source of the
long-chain conclusion.

## What cannot be claimed

- exact equality or a mathematical proof of the FT;
- symmetry for every `k` in `[0,1]`;
- uniform accuracy as `n` tends to infinity;
- that action current itself is entropy production;
- that the medium-entropy direct-sampling slopes at `t<=200` already reached
  their asymptotic limit;
- a universal statement outside the simulated bath parameters and numerical
  dynamics.

## Traceability

- final numerical table: `FINAL_SUMMARY.csv`;
- prospective protocol and matrix: `PROTOCOL.md`, `RUN_MATRIX.csv`;
- raw Phase-II output: `raw/`;
- accepted analysis: `analysis/final/`;
- analysis-grid erratum: `ANALYSIS_ERRATUM.md`;
- full audit: `VALIDATION_REPORT.md`;
- mentor direct-sampling audit:
  `../entropy_ft_2026-08-26/production/final_v4/`;
- `n=2` total-entropy audit:
  `../entropy_ft_scgf_2026-08-27/total_entropy_n2_short/parametric_analysis/`;
- independent discrete path-ratio audit:
  `../discrete_path_ft_2026-08-28/production_v2/`.
