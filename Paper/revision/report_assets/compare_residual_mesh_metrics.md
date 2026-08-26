# LTE residual mesh metrics

Generated: `2026-06-21T07:05:01.565639+00:00`

These metrics reproduce the plotting convention of
`/Users/jayleenjiang/Documents/MATLAB/lte/compare_residual.m` for
`report_assets/compare_residual_mesh.pdf`.  They are descriptive
mesh-slice diagnostics and are not replacements for the fixed
weighted-core LTE estimators in `source_trace_metrics.json`.

## Method

- Histogram grid: `NB=80`.
- MATLAB theta index: `20`.
- Fit mask: `Q>50` and `P>50` over the full histogram.
- Surface display mask: `Q>=50` and `P>0` on the displayed slice.
- Reported norms are unweighted over displayed mesh bins.

## Summary

| case | pair | slope x | displayed bins | slice RMS | slice mean abs | core RMS | residual range |
|---|---:|---:|---:|---:|---:|---:|---|
| `n15` | 4 | 0.795614 | 4809 | 0.253803 | 0.226076 | 0.248066 | [-0.822048, 0.595117] |
| `n25` | 6 | 0.780028 | 4040 | 0.199409 | 0.169806 | 0.193853 | [-0.729536, 0.944319] |
| `n50` | 12 | 0.840796 | 2891 | 0.137960 | 0.089406 | 0.116037 | [-0.560516, 1.292206] |
