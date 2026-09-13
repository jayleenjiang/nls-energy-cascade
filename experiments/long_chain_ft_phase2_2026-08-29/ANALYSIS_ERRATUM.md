# Analysis-script erratum

The prospectively frozen simulation matrix, estimands, final horizons,
thresholds, and conditional stopping rule are unchanged.

The first post-run invocation of `analyze_long_chain_ft.py` stopped before
producing an accepted analysis because the wrapper requested horizons 30 and
50 for data recorded every two time units.  The analyzer also requires an
exact late-half row, but times 15 and 25 are absent from those timeseries.
This was a discrete-indexing error in the wrapper, not a simulation failure.

The failed partial output is retained under `analysis/n30_population/`.
Accepted analyses are written under `analysis/final/` and use horizons
20, 40, and 60 for `n=20,30`.  Their exact half-horizons 10, 20, and 30 are
recorded.  The frozen primary decision remains at horizon 60, and every
pre-registered numerical and convergence threshold is unchanged.  The
`n=40` horizons were already compatible with the recorded grid and are not
changed.
