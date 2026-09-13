# Validation report

- Protocol frozen before analysis: commit `0af6652`.
- Analysis source frozen before output: commit `70cbab9`.
- No simulation was launched; four existing controlled trajectory ensembles
  were reused.
- Inputs: 4 cases x 64 streams; all 256 files have expected size.
- Pair quadrature: 2048 fixed origins per stream, 131072 pairs per case.
- Dictionaries: nominal E1/E2/E3 = 109/269/519; retained dimensions are listed
  in the PDF and raw Gram diagnostics.
- Lags: short 0.02/0.05/0.10; long 0.25/0.50/1.00.
- Nyquist exclusion: `abs(Im lambda) >= 0.8*pi/tau` marked aliased.
- Gram cutoff: primary `1e-10`; sensitivities `1e-8`, `1e-12`.
- Constant-mode gate: all errors <= `1.28e-12`, below `1e-8` threshold.
- Bootstrap: 500 whole-trajectory resamples per case for the oscillatory target
  and 500 per tracked real candidate.
- Oscillatory gate: 3 PASS, 1 FAIL (driven `dt=1e-3` dictionary spread).
- Leading visible real candidates: all local convergence gates pass, but the
  result moves from the old dictionary in 3/4 interval comparisons.
- Integrity audit: `AUDIT.json` status PASS.

The statistical and provenance audit passes. The scientific completeness gate
does not.
