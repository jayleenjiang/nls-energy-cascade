# Entropy SCGF and rare-event FT study

This experiment follows the audited direct two-tail study in
`../entropy_ft_2026-08-26/`.  Its first gate is a zero-new-simulation pilot:
estimate the finite-time medium-entropy SCGF from the existing million-sample
blocks and quantify where exponential reweighting loses effective support.

The target observable is

```text
psi_t(k) = log E[exp(-k Sigma_medium)] / t.
```

The asymptotic Gallavotti--Cohen diagnostic is

```text
psi(k) = psi(1-k).
```

At finite time the medium entropy lacks the system-entropy endpoint term, so
the symmetry is not expected to be exact.  Direct estimates are accepted only
where sample-weight and independent-stream ESS gates pass.  The completed
cloning pilot validates the low-tilt point at `n=10`, but the complementary
high-tilt point and the longer-time extension fail ancestry-support gates.
Consequently no `n=20,30,40` cloning production was launched.

## Direct-pilot command

```sh
python3 flux/analyze_entropy_scgf.py \
  experiments/entropy_ft_2026-08-26/production/n10_blocks.csv \
  experiments/entropy_ft_2026-08-26/production/n20_blocks.csv \
  experiments/entropy_ft_2026-08-26/production/n30_blocks.csv \
  experiments/entropy_ft_2026-08-26/production/n40_blocks.csv \
  --output-dir experiments/entropy_ft_scgf_2026-08-27/direct_pilot

python3 flux/audit_entropy_scgf.py \
  experiments/entropy_ft_scgf_2026-08-27/direct_pilot
```

Raw files are read-only inputs and are not copied or modified.
Their SHA-256 values, together with the endpoint-control raw files and all
analysis programs, are recorded in `SOURCE_MANIFEST.sha256`.  The large raw
CSV files remain local/ignored; the audited summaries, diagnostics, figures,
and hashes are versioned.

## Cloning build and checks

```sh
clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  flux/NLS_entropy_cloning.cpp -o entropy_cloning

./entropy_cloning selftest
```

The accepted low-tilt population check uses four independent runs at
`n=10`, `t=20`, `k=0.3`, `N_c=1024`, `Delta=0.5`, and `dt=5e-4`.
Aggregate and audit all available cloning summaries with

```sh
python3 flux/analyze_entropy_cloning.py \
  'experiments/entropy_ft_scgf_2026-08-27/**/*summary*.csv' \
  --direct-scgf \
  experiments/entropy_ft_scgf_2026-08-27/direct_pilot/scgf_points.csv \
  --output-dir experiments/entropy_ft_scgf_2026-08-27/aggregate

python3 flux/audit_entropy_cloning.py \
  experiments/entropy_ft_scgf_2026-08-27 \
  experiments/entropy_ft_2026-08-26/production/n10_summary.csv \
  --output experiments/entropy_ft_scgf_2026-08-27/aggregate/audit.md
```

`PILOT_STATUS.md` states the allowed conclusion.  In particular, successful
low-tilt validation does not open the FT gate while `k>=0.5` lacks ancestry
support.

## Separate n=2 total-entropy endpoint pilot

The endpoint mode writes the reduced start/end state and thermodynamic block
observables without modifying the existing production sampler:

```sh
./entropy_cloning endpoints T_left T_right 2 streams burnin block_time \
  blocks_per_stream dt seed threads out_prefix
```

An equal-temperature data set selects and validates the density smoothing;
the driven data are held out from that choice.

```sh
python3 flux/analyze_total_entropy_n2.py \
  equilibrium_blocks.csv driven_blocks.csv \
  --equilibrium-temperature 6 \
  --output-dir total_entropy_n2/analysis

python3 flux/audit_total_entropy_n2.py total_entropy_n2/analysis \
  --output total_entropy_n2/analysis/audit.md
```

See `TOTAL_ENTROPY_CONTROL_DESIGN.md` and the predeclared gates before
interpreting this calculation.  The completed control passes pointwise
equal-temperature density-ratio checks but is blocked by the learned
equilibrium integral-FT test and a driven exponential-weight ESS of only
`1.63`.  The accepted overall statement is in `PILOT_STATUS.md`.
