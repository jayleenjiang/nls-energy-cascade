# Canonical action-current experiment

`NLS_flux_canonical.cpp` is the publication-oriented replacement for the
June 18 flux program.  The old source and data are retained for provenance but
must not be used as evidence for the manuscript because their boundary-noise
amplitudes omit the required factor `sqrt(2)`.

The canonical program:

- implements only the Gibbs-preserving boundary bath;
- labels the measured quantity as the **action current**, not heat/energy;
- uses fixed-step Euler--Maruyama in double precision;
- uses standard trigonometric and logarithmic functions;
- assigns a deterministic independent RNG stream to every trajectory;
- writes raw finite-time trajectory averages and four equal current blocks,
  action profiles, burn-in
  diagnostics, complete run metadata, and positivity-projection diagnostics;
- reports the standard error `s/sqrt(N)` and a normal 95% confidence interval.

Build on Apple Silicon:

```sh
clang++ -O3 -mcpu=native -std=c++17 \
  -Xpreprocessor -fopenmp \
  -I/opt/homebrew/include/eigen3 \
  -I/opt/homebrew/opt/libomp/include \
  -L/opt/homebrew/opt/libomp/lib -lomp \
  NLS_flux_canonical.cpp -o flux_canonical
```

Run:

```sh
./flux_canonical T1 Tn n batches burnin measure dt seed threads out_prefix [bond]
```

One batch contains 16 trajectories.  The optional bond is the right endpoint
`j` of `(j-1,j)` and defaults to `n/2`.

For each run, the following files are written:

- `<prefix>_samples.csv`
- `<prefix>_summary.csv`
- `<prefix>_profile.csv`
- `<prefix>_burnin.csv`

The sample file contains the full-window current, the two half-window currents,
and four equal block currents.  A run with `measure=200` therefore supports
direct comparisons at averaging windows 50, 100, and 200 without rerunning the
simulation.

Analyze a set of runs:

```sh
python3 analyze_canonical_flux.py results/*_summary.csv \
  --primary-dt 0.0005 \
  --bootstrap 10000 \
  --seed 20260619 \
  --output-prefix results/flux
```

The bootstrap interval quantifies Monte Carlo uncertainty conditional on the
chosen chain lengths, timestep, burn-in, averaging window, and power-law model.
It does not include finite-size or discretization-model uncertainty; those are
reported separately through timestep and window sensitivity.
