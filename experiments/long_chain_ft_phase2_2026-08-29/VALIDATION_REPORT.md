# Validation report and material passport

## Material passport

### Scientific question

The direct-sampling study asks how both current tails depend on threshold and
averaging time.  The thermodynamic FT question asks whether total entropy at
finite time, or its long-time entropy-current SCGF, has the required
fluctuation symmetry.  Action-current statistics are retained as a separate
observable and are not substituted for entropy production.

### Dynamics and observables

- Projection-free Cartesian variables `c_j=x_j+i y_j`; no action floor.
- Baths `(T_L,T_R)=(10,2)`, coupling `gamma=0.1`.
- Bath heat is accumulated from the bath increment, not from total energy
  change; medium entropy uses
  `Sigma_m=-Q_L/T_L-Q_R/T_R`.
- The long-chain asymptotic observable is the right-bath entropy current in
  the gauge `Sigma_R=-0.4 Q_R`.

### Code provenance

- Simulator: `flux/NLS_entropy_cloning.cpp`, pre-run SHA-256
  `bc3c27bf62a45aa879c4a8fd3e4d70fe7bdf1e3bcb0d06c0318826e61e615b6d`.
- Analysis engines: `flux/analyze_long_chain_ft.py` and
  `flux/analyze_long_chain_ft_controls.py`; their pre-run hashes remain
  unchanged.
- Frozen parent revision: `d37f5336ae43b1a1326b281a387f3a346993a000`.
- Prospective protocol commit: `32c86d3`; analysis-grid erratum commit:
  `a0bee1f`.
- The controlled Gaussian-kernel and cloning self-tests passed before
  production.

The pre-run manifest now reports one expected mismatch: the wrapper
`analyze_phase2.sh` was corrected after its first invocation requested
late-half rows at unrecorded odd times.  `ANALYSIS_ERRATUM.md` preserves the
failed partial output and explains the correction.  The executable source,
analysis engines, final horizons, estimands, thresholds, and simulation data
were not changed.

### Data passport

| material | size | purpose | accepted audit |
|:---|:---|:---|:---|
| mentor direct blocks | 1,000,064 blocks for each `n=10,20,30,40` | two tails and `t=20,...,200` dependence | 75,617 independent checks, 0 errors |
| finite-time `n=2` NESS blocks | 1,048,576 per condition | total-entropy detailed and integral FT | all parametric gates pass |
| discrete path-ratio data | 1,000,000 forward/reverse paths at accepted step | independent Crooks/IFT control | all primary gates pass |
| Phase-II cloning output | 56 summaries, 56 timeseries, 56 logs | missing long-chain convergence controls | 0 malformed/nonfinite/incomplete files |

Phase-II uses four new independent seeds for each tilt member in every cell.
No failed cell was followed by data-dependent retuning.

## Recomputed Phase-II audit

- Raw completion: 56/56 summaries, 56/56 timeseries, 56/56 logs.
- Nonfinite numeric entries: zero.
- Incomplete timing logs: zero.
- Midpoint failures: zero.
- Every final pair passes the minimum selection-weight ESS requirement.
- Final symmetry intervals contain zero and absolute residuals are at most
  0.01 at `n=10,20,30,40`; late-half conclusions agree.
- Observation-time convergence passes for every accepted baseline and
  control group.
- Population convergence passes for both `n=30` members and their paired
  residual from `N_c=2048` to 4096.
- At each of `n=20,30,40`, all four member controls and both paired timestep
  and selection controls pass.  The earlier `n=10` control suite also passes.

Exact machine-readable values are in `FINAL_SUMMARY.csv` and the generated
CSVs under `analysis/final/`.

## Direct-sampling cross-check

The four raw direct files retain their accepted SHA-256 hashes.  The raw
audit confirms exactly 1,000,064 ordered rows per chain length, no midpoint
failures, and energy-balance RMS rates from `9.71e-6` (`n=10`) to `3.05e-6`
(`n=40`).  The independent analysis audit recomputed all reported aggregates,
tail counts, probability ratios, fits, stationarity statistics, and
heat--action coupling metrics with 75,617 checks and zero errors.

## Eleven-fallacy scan

1. **Non-rejection as proof:** avoided; passing is described as numerical
   consistency, not exact equality.
2. **Data-dependent tilt selection:** avoided; complementary pairs and gates
   were frozen before Phase-II output.
3. **Optional stopping/retuning:** avoided; the `n=30` branch and its stopping
   rule were prospective.
4. **Pseudoreplication:** avoided in the cloning comparison by four
   independent seeds per member; direct blocks retain stream identifiers.
5. **Unsupported-tail smoothing:** avoided; plus-four and Gaussian fits are
   descriptive and never manufacture a resolved FT tail.
6. **Finite-time/asymptotic conflation:** avoided; direct `t<=200` and
   long-time SCGF statements are separated.
7. **Medium versus total entropy:** avoided; the endpoint system entropy is
   included only in the finite-time `n=2` total-entropy result.
8. **Action current as entropy:** avoided; imperfect heat--action coupling is
   reported explicitly.
9. **Discretization blindness:** avoided by timestep-halving controls and the
   independent discrete-kernel convergence audit.
10. **Absolute-only effect reporting:** mitigated by reporting the normalized
    SCGF residual; the larger relative uncertainty at `n=40` is retained.
11. **Finite-size extrapolation:** avoided; no infinite-chain or all-tilt
    conclusion is drawn.

## Reproducibility status

**Reproducible from archived code and seeds.**  `ENVIRONMENT.md` records the
Apple M4/Clang/OpenMP/Python environment.  `build_and_selftest.sh`,
`run_phase2.sh`, and `analyze_phase2.sh` reproduce the Phase-II workflow.
Raw and generated Phase-II files are small enough to remain in the Git
history.  The much larger million-block direct and `n=2` raw files remain in
their documented local experiment directories; their accepted hashes and
audits provide integrity anchors.

## Final scope assessment

The evidence is sufficient for a carefully scoped numerical-paper claim that
the resolved entropy-current SCGF pairs are consistent with GC symmetry at
all four tested finite chain lengths.  A stronger claim would require a
predeclared grid of additional tilts and, for an infinite-chain statement, a
separate finite-size convergence study of the SCGF symmetry residual itself.
