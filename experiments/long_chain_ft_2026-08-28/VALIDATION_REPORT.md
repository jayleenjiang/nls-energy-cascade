## Material Passport

- Origin Skill: experiment-agent
- Origin Mode: validate
- Origin Date: 2026-08-29
- Verification Status: VERIFIED
- Version Label: long_chain_ft_validation_v1

## Validation Report

- **Source**: `long_chain_ft_2026-08-28`
- **Overall Confidence**: CAUTION

`VERIFIED` means that the code path, deterministic recovery, analysis
regeneration, and frozen audit logic were successfully rerun in the recorded
environment.  `CAUTION` describes the scientific scope: the broad all-chain FT
hypothesis did not pass every pre-specified numerical gate.

### Statistical Findings

| Metric | Test | Value | Effect size / tolerance | Confidence |
|:---|:---|:---|:---|:---|
| `n=10` GC residual | independent-run Welch interval, frozen paired/time/population gates | `D=-0.000075`, 95% CI `[-0.009953,0.009802]` | `|D|<0.01`; timestep and supported selection controls pass | SOLID within tested pair |
| `n=20` GC residual | same core audit | `D=-0.003389`, 95% CI `[-0.012081,0.005302]` | core gates pass; strongest `(0.3,0.7)` pair unresolved | CAUTION |
| `n=30` GC residual | same core audit at `N_c=2048` | `D=-0.003574`, 95% CI `[-0.008235,0.001087]` | pair/time/support pass; `k=0.4` population change is `-0.003470 = 2.36` combined SE | CAUTION / unresolved |
| `n=40` GC residual | independent remedial series at `t=120`, `N_c=1024` | `D=0.004498`, 95% CI `[-0.002343,0.011340]` | core and timestep gates pass; selection control fails | CAUTION / unresolved |
| `n=40` selection sensitivity | supported `selection_time 2 -> 1`, `N_c=512` | max member change `0.005727`; late residual change `0.010309` | member is `2.84` combined SE; late change exceeds `0.01` | RED_FLAG for full-control claim |
| gauge identity | direct million-block samples | RMS `6.1e-6`--`1.9e-5` | numerical first-law/gauge discrepancy | SOLID |
| local path-ratio control | arbitrary-chain forward/reverse transient tests | supported controls pass; independent `n=30`, `N=5e5` Crooks slope `0.9941` | target slope 1 | SOLID for tested transient kernel |

The intervals quantify variation among four independent simulation seeds.  The
absolute `0.01` and two-combined-SE rules are predeclared numerical acceptance
gates, not p-value significance tests and not a multiple-comparison-corrected
proof.

### Warnings

| Type | Detail | Affected |
|:---|:---|:---|
| finite numerical scope | Only one resolved nontrivial pair is certified at each reported length; `(0.3,0.7)` loses support for `n>=20`. | any all-`k` claim |
| small number of independent runs | Every aggregate uses four independent seeds; Student/Welch intervals are therefore wide and sensitive to one run. | all SCGF intervals |
| population convergence | The final `n=30`, `k=0.4` member fails the unchanged two-combined-SE population gate. | `n=30` final status |
| estimator sensitivity | The supported `n=40` selection control fails; the earlier selection-4 control also loses ESS support. | `n=40` and all-chain status |
| endpoint-gauge condition | Right-current and medium entropy differ by an endpoint energy term; equal asymptotic SCGFs require controlled endpoint exponential moments. | theorem-level interpretation |
| finite-size/asymptotic boundary | No analytic ergodicity/domain proof and no `n -> infinity` result is supplied. | “proof” or thermodynamic-limit wording |

### Fallacy Scan

- **Coverage**: 11/11 fallacy types checked

| Fallacy | Severity | Detail | Recommendation |
|:---|:---|:---|:---|
| Simpson's paradox | NOTE | No subgroup aggregation is used; each `n,k,N_c,t` group is audited separately. | N/A |
| Ecological fallacy | NOTE | No individual-level inference is made from group aggregates. | Keep claims at the simulated chain/group level. |
| Berkson's paradox | NOTE | No outcome-conditioned sampling frame is used; rare-event control is likelihood-corrected. | Retain exact likelihood-ratio self-tests. |
| Collider bias | NOTE | No regression adjustment or conditioned collider enters the SCGF comparison. | N/A |
| Base-rate neglect | NOTE | This is not a diagnostic-classification analysis. | N/A |
| Regression to the mean | NOTE | Production groups are not selected by extreme outcome values; independent remedial seeds are not pooled with failures. | Continue reporting original and remedial series separately. |
| Survivorship bias | NOTE | Failed strong-pair, horizon, population, ESS, and selection diagnostics are retained rather than deleted. | Preserve raw failures and the claim ledger. |
| Look-elsewhere effect | CAUTION | Multiple tilts, horizons, controls, and pilot settings were explored.  Mitigation is predeclared secondary pairs, frozen gates, and explicit failed-run retention. | Do not use unregistered pairs as confirmatory evidence. |
| Garden of forking paths | CAUTION | Remedial runs followed observed failures.  Each amendment was frozen before new seeds and prohibited further chasing after a repeated failure. | Present amendments and stopping rules with the result. |
| Correlation does not imply causation | NOTE | The claim concerns a dynamical symmetry, not an observational causal effect. | Avoid unrelated causal wording. |
| Reverse causality | NOTE | No directional observational association is interpreted. | N/A |

### Reproducibility

- **Method**: deterministic same-seed recovery plus stochastic independent-seed
  regeneration in the recorded Apple M4/Clang/OpenMP/Python environment
- **Verdict**: REPRODUCIBLE for the scoped code-and-audit pipeline

| Metric | Original | Re-run | Diff | Status |
|:---|:---|:---|:---|:---|
| recovered `n=40`, `N_c=512` aggregate (6 rows x 6 fields) | archived aggregate | same-seed recovery | maximum absolute difference `0` | MATCH |
| zero-control change of measure | original finite-step estimator | controlled estimator at zero shift | machine precision | MATCH |
| Gaussian cloning self-test | exact symmetry target | release and sanitizer builds | pass | MATCH |
| final core audit | four frozen rows | regenerated from raw summaries | 3/4 core rows pass | MATCH |
| endpoint numerical audit | four frozen controls | regenerated from raw summaries | 3/4 controls pass | MATCH |

The reproducibility verdict does not override failed scientific gates: a
reproducible unresolved result remains unresolved.
