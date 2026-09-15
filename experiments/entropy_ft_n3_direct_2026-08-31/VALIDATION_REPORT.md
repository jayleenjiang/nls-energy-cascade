# Validation report

## Reproducibility status

`VERIFIED` for data integrity and deterministic analysis rerun; `NOT
ESTABLISHED` for the physical finite-time FT because direct negative-tail
support fails and no admissible `(10,2)` stationary density is available.

## Evidence map

| item | evidence | status |
|---|---|---|
| Frozen scientific parameters and gates | `PROTOCOL.md`, Git commit `9fd7d75` | verified before production |
| Source extension | `flux/NLS_entropy_ft.cpp`, SHA-256 `9ae5835e...6d4eaf2` | self-test and legacy-mode comparison pass |
| Production command and binary | `provenance/production_manifest.txt` | exact command, binary hash, seed recorded |
| Raw row integrity | `analysis/analysis_audit.json` | 1,000,064 x 19, all finite and ordered |
| Heat/entropy identity | same audit | pass within scale-aware roundoff |
| First-law residual | `analysis/first_law_residuals.csv` | fully reported, no row deletion |
| Negative-tail counts | `analysis/negative_tail_counts.csv` | 41 at t=20; zero at t>=40 |
| Symmetric raw bins | `analysis/symmetric_bin_counts.csv` | zero qualifying pairs |
| Medium exponential average | `analysis/medium_entropy_ift.csv` | numerically unresolved; ESS about 1--5 |
| Section-4 NN compatibility | `MODEL_COMPATIBILITY_AUDIT.md` | fails parameter and precision requirements |
| Initial analysis failure | `analysis_failed_v1/`, `ANALYSIS_ERRATUM.md` | preserved; analysis-only repair documented |

## Statistical-fallacy scan

1. **Outcome-dependent fit selection:** avoided; no slope is fitted because the
   frozen raw-count gate fails.
2. **Treating zero observed events as zero probability:** avoided; the report
   states raw zeros and does not infer an exact zero tail.
3. **Confusing medium with total entropy:** avoided explicitly.
4. **Calling absence of support an FT violation:** avoided explicitly.
5. **Ignoring temporal dependence:** aggregation stays within stream and CIs
   resample independent streams.
6. **Pseudo-replication:** block counts and stream count are both disclosed.
7. **Tail-dominated exponential means:** ESS and maximum weight share reported.
8. **Post-hoc bin/window tuning:** avoided; FD bins and gates were frozen.
9. **Numerical-integrator error omission:** first-law distribution and
   midpoint failures reported.
10. **Model-domain extrapolation:** avoided; the `(2,8)` NN is not applied to
    `(10,2)` endpoints.
11. **Positive-result bias:** the negative support result is retained and no
    cloning output is substituted.

## Strict claim boundary

The experiment establishes only that one million ordinary `t=20` blocks are
insufficient to resolve the `n=3` two-sided medium-entropy distribution under
the frozen gate, and that longer direct aggregations contain no observed
negative events.  It does not establish an exact FT, an FT violation, or a
total-entropy result.
