# Long-chain Gallavotti--Cohen verification

This experiment targets the long-time SCGF symmetry for the boundary-driven
Cartesian NLS chain at `n=10,20,30,40`:

```text
psi_n(k) = lim_(t->infinity) t^(-1) log E_NESS exp[-k Sigma_R(t)],
psi_n(k) = psi_n(1-k).
```

The primary additive observable is the right-bath entropy current

```text
Sigma_R = -(beta_R-beta_L) Q_R.
```

Here `Q_R` is heat in the Cartesian physical energy `E=H/2`.  Consequently
`rho_eq proportional exp(-E/T)=exp(-H/(2T))`, exactly matching the manuscript's
Gibbs normalization.  Equivalently, if heat is expressed in the manuscript
Hamiltonian units `Q_R^H=2 Q_R`, the inverse-temperature difference is `0.2`
rather than `0.4`; the entropy current is unchanged:
`-0.4 Q_R=-0.2 Q_R^H`.

With heat positive into the system,

```text
Sigma_R - Sigma_medium = beta_L Delta E + O(first-law residual).
```

Thus both observables have the same long-time SCGF when the endpoint-energy
exponential moments are controlled.  `Sigma_R` is used for rare-event
sampling because it removes most of the cold-bath endpoint-energy term that
destroyed direct and medium-entropy cloning support near `k=1`.

The controlled cloning sampler uses a modified bath drift only as an
importance proposal.  Every finite-step proposal is corrected by the exact
Gaussian transition-density likelihood ratio, so no fluctuation symmetry is
assumed in the estimator.

## Authoritative programs

- `flux/NLS_entropy_cloning.cpp`: `controlled` mode implements exact finite-step
  controlled importance sampling plus fixed-population cloning.
- `flux/analyze_entropy_current_scgf.py`: direct SCGF, gauge-identity, support,
  and finite-time symmetry analysis.
- `flux/analyze_long_chain_ft.py`: independent-run aggregation and acceptance
  audit (added after the pilot design is frozen).
- `flux/analyze_long_chain_ft_controls.py`: fixed-timestep and
  selection-interval comparison at a common population and horizon.
- `flux/summarize_long_chain_ft.py`: compiles only the rows frozen in
  `FINAL_RESULT_SPEC.csv` into the final evidence table and figure; it exits
  unless every required row exists.
- `flux/summarize_long_chain_ft_controls.py`: compiles the endpoint timestep
  and selection-interval comparisons into the report table.

## Status

The existing million-block direct data are retained as a baseline.  New
controlled-cloning results are accepted only over the chain lengths, horizons,
and paired `k` values that pass the gates in `PREDECLARED_GATES.md`.

The final frozen audit is intentionally mixed.  The complete core and
numerical-control suite passes at `n=10`.  The core paired-SCGF audit passes at
`n=20` and `n=40`, but the endpoint selection audit fails at `n=40`, so no
fully controlled all-long-chain claim is made.  At `n=30`, the paired residual
is consistent with zero but an individual population-convergence gate fails.
See `FINAL_VERDICT.md` for the claim boundary and `VALIDATION_REPORT.md` for
the reproducibility and 11-item fallacy audit.

The first four-seed `n=10`, `N_c=1024` diagnostic used
`control_scale=0.25`.  It was rejected before extension to other chain lengths:
one `k=0.7` run crossed the predeclared weight-ESS floor and the cumulative
paired residual had not converged by `t=40`.  A support-only refinement scan
(weight ESS, root-weight ESS, and genealogy; the GC residual was not used for
selection) compared scales `0.25`, `0.50`, and `0.75`.  Scale `0.50` was frozen
for the replacement production series because it removed the low-weight-ESS
failure without the severe state/genealogy collapse seen at `0.75`.  See
`PILOT_DECISIONS.md` for the numerical audit trail.

A subsequent support-only selection-interval audit rejected the original
`selection_time=0.5`: repeated neutral resampling caused avoidable genealogical
collapse at long horizons.  The replacement production series therefore uses
`selection_time=2.0`, which recovered the independently simulated direct
`t=20` SCGF within uncertainty while preserving substantially more roots.  No
GC residual was used to select this interval.

The strongest pair `(0.3,0.7)` is resolved for `n=10`.  For `n>=20`, its
high-tilt member is support/finite-time limited, so it is retained as a failed
diagnostic rather than reported as a positive result.  The longer chains are
tested on the already predeclared secondary pair `(0.4,0.6)`.  Any final claim
must state this narrower resolved tilt interval explicitly.

The first frozen cross-chain summary subsequently exposed failed individual
`N_c=512 -> 1024` population-member gates at `n=30`, even though the paired
residual gate passed.  The failure is retained.  Before any further data were
generated, `REMEDIAL_AMENDMENT_2026-08-29.md` froze one independent
`N_c=2048` extension and its unchanged acceptance criteria.

## Reproduction order

Build the release sampler with the recorded Apple-Clang/OpenMP environment:

```bash
clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  flux/NLS_entropy_cloning.cpp -o flux/entropy_cloning_v2
flux/entropy_cloning_v2 selftest
```

The frozen execution order is:

```bash
bash experiments/long_chain_ft_2026-08-28/run_production.sh all
bash experiments/long_chain_ft_2026-08-28/run_n30_population_extension.sh
bash experiments/long_chain_ft_2026-08-28/run_n40_remedial_extension.sh all
bash experiments/long_chain_ft_2026-08-28/run_numerical_controls.sh all
bash experiments/long_chain_ft_2026-08-28/run_numerical_controls_n10_remedial.sh all
bash experiments/long_chain_ft_2026-08-28/run_numerical_controls_n40.sh all
bash experiments/long_chain_ft_2026-08-28/run_selection_control_n40_remedial.sh
bash experiments/long_chain_ft_2026-08-28/analyze_all.sh
bash experiments/long_chain_ft_2026-08-28/analyze_controls.sh
```

The final pre-specified tables and figure are compiled only after the audit
directories exist:

```bash
python3 flux/summarize_long_chain_ft.py \
  experiments/long_chain_ft_2026-08-28 \
  --output-dir experiments/long_chain_ft_2026-08-28/final_summary
python3 flux/summarize_long_chain_ft_controls.py \
  experiments/long_chain_ft_2026-08-28 \
  --output-dir experiments/long_chain_ft_2026-08-28/final_controls
```
