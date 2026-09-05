# Validation log

## Source and binary checks

- Release build: Apple clang 17, C++17, `-O3 -mcpu=native`, OpenMP.
- `entropy_cloning_v2 selftest`: controlled zero-shift kernel identity PASS;
  Gaussian cloning symmetry self-test PASS.
- Address/undefined-behavior sanitizer build: self-tests and a controlled
  `n=4` smoke run PASS with no sanitizer diagnostics.
- `discrete_path_ft selftest`: midpoint reversibility, Hamiltonian
  time-reversal, conjugate-path recovery, and bath-kernel ratio PASS.
- Address/undefined-behavior sanitizer build of the path-ratio code: self-test
  and `n=10` smoke run PASS with no sanitizer diagnostics.
- Python analysis scripts compile with `python3 -m py_compile`.
- Energy normalization audit PASS: the measured heat is the increment of
  `E=H/2`; therefore the code's `-Q_E/T` entropy is the same convention as the
  manuscript Gibbs density `exp[-H/(2T)]`.

## Independent numerical controls

- Burn-in `B=500` is inherited from the separate profile/current relaxation
  study for `n<=40`; the original no-burn-in evidence remains in the burn-in
  experiment module.
- Million-block direct estimates calibrate the controlled SCGF at `t=20`
  wherever direct exponential averages retain support.
- The right-current/medium-entropy gauge identity has RMS
  `6.1e-6`--`1.9e-5` over `n=10,20,30,40`.
- At low resolved tilts, the directly sampled finite-time gauge-SCGF
  difference decreases with observation time.  It is retained as a finite-time
  endpoint effect, not forced to zero.
- Arbitrary-chain one-step path-ratio controls pass the transient Crooks/IFT
  gates for `n=10,20,40`; a separately seeded `N=5e5` replication passes for
  `n=30`.  Failed longer-window/low-sample attempts remain on disk.
- Two `n=40`, `N_c=512` run-3 files that macOS/iCloud had evicted were
  materialized again.  A separate deterministic recovery with the same seeds
  reproduced 6 horizon/tilt rows and 6 diagnostic fields per row exactly
  (`max_abs_diff=0`).  No original output was overwritten.
- The first `n=10`, `N_c=512` numerical-control series is retained as a failed
  diagnostic: the baseline high-tilt member missed the fixed ESS support floor,
  and selection interval 4 had still lower support and no sample at the
  required half horizon.  A replacement matrix was frozen before new data were
  generated.  At `N_c=1024`, timestep halving (`5e-4 -> 2.5e-4`) and selection
  interval shortening (`2 -> 1`) both pass all unchanged support, member,
  full-window-residual, and late-half-residual comparison gates.

## Long-chain production status

- `n=10`, pair `(0.3,0.7)`: paired, support, time, and population gates pass
  through `t=60`.
- `n=20`, pair `(0.4,0.6)`: paired, support, time, and population gates pass
  through `t=60`; `(0.3,0.7)` remains unresolved.
- `n=30`, pair `(0.4,0.6)`: the `N_c=2048` paired residual is `-0.003574`
  (95% interval `[-0.008235,0.001087]`) and its support, paired, and time gates
  pass.  However, the frozen 1024-to-2048 `k=0.4` individual-member change is
  `-0.003470`, about `2.36` combined SE, so the unchanged population gate
  fails.  Per the pre-run amendment no additional population is selected;
  `n=30` remains unresolved.  The stronger `(0.3,0.7)` pair is also unresolved.
- `n=40`, pair `(0.4,0.6)`: the first `N_c=1024` series passes support,
  full-window CI, population, and `t=60` to `t=80` time comparisons, but fails
  the unchanged late-half absolute gate (`0.01527` at `t=60`, `0.01013` at
  `t=80`).  It is retained as a failed diagnostic.  The independently seeded
  `t=120`, `N_c=512/1024` remedial series passes the paired, late-half, time,
  and population gates.  At `N_c=1024`, the final residual is `0.004498`
  (95% interval `[-0.002343,0.011340]`) and the late-half residual is
  `0.007935` (95% interval `[-0.005986,0.021855]`).
- Endpoint numerical controls: all `n=10`, `N_c=1024` timestep and supported
  selection comparisons pass.  At `n=40`, timestep halving passes, while the
  original selection-4 series fails support.  The pre-specified supported
  selection-1 replacement also fails: its `k=0.4` member changes by `0.005727`
  (`2.84` combined SE), and the late-half paired residual changes by `0.010309`,
  just beyond the unchanged `0.01` absolute gate.  No further interval is
  selected.  The `n=40` core GC result is therefore not promoted to a fully
  controlled positive claim.

## Remaining before final claim

- regenerate final summary tables/figures;
- compile and visually inspect the report;
- run source/data manifests, `git diff --check`, and repository sync.
