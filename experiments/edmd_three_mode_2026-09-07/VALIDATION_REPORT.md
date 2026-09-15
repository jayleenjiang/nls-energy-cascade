## Material Passport

- Material ID: `NLS-EDMD-N3-2026-09-07`
- Type: computational experiment / validation report
- Status: `VERIFIED` for reproducibility and gate application; scientific
  outcome is a partial positive result with a failed oscillatory known-answer
  check.
- Inputs: four completed controlled-integrator trajectory ensembles, 64
  streams each, 100001 snapshots per stream.
- Analysis: frozen EDMD protocol, 131072 lagged pairs per case, three nested
  dictionaries, four lags, two timesteps, 500 whole-stream bootstraps.

## Result

The constant mode, Gram conditioning, cutoff sensitivity, and output-integrity
checks pass.  In the driven system only one nontrivial eigenvalue is reportable
at both timesteps, with rates -0.940174 and -0.963191.  The next candidate
rates cannot be separated by their bootstrap intervals.  The 5.3 oscillation
fails lag convergence.  The EDMD study therefore does not establish several
distinct driven slow modes; it supports one leading slow mode plus unresolved
faster modal contamination.

At equilibrium, rates near -0.95 and -2.85 are mutually reproduced at both
timesteps and pass the frozen gates.

## Reproducibility verification

- Raw trajectory data were opened read-only; no simulation was run.
- The 5D reduced state was reconstructed algebraically from saved columns.
- Every expected spectrum, bootstrap, Gram, matrix, slice, and figure artifact
  was found.  The automated audit reports 42/42 checks passed.
- All 20 D3 matrices have shape 500 x 500 and are stored as compressed CSV.
- Every reference mode has exactly 500 bootstrap rows; matches are retained in
  the raw bootstrap tables.
- Full artifact hashes are in `provenance/analysis_outputs.sha256`.

## Statistical and inferential cautions

1. Time origins within a stream are correlated; uncertainty is therefore
   resampled over whole independent streams, not individual origins.
2. EDMD modes are non-normal.  The reported observable weights are signed modal
   covariance contributions, not probabilities.
3. The oscillatory known-answer failure limits the claim: the dictionary is not
   adequate for the complete accessible spectrum.
4. The faster driven candidates overlap statistically.  Reporting them as two
   eigenvalues would be over-resolution.
5. The finite dictionary approximates Koopman modes; no exact generator theorem
   is inferred from this numerical analysis.

## Fallacy scan

- No p-value dichotomisation: checked.
- No CI-as-proof interpretation: checked.
- No post-hoc dictionary, lag, threshold, or mode-window tuning: checked.
- No pseudo-replication of time origins: whole-stream bootstrap used.
- No causal claim from modal association: checked.
- No extrapolation beyond n=3 or the two bath settings: checked.
- No selective omission of failed modes: all spectra and failed gates retained.
- No multiple-comparison claim based on point estimates alone: neighbour-CI
  gate applied.
- No conflation of numerical stability with physical truth: oscillation control
  failure is explicit.
- No equivalence claim from failure to reject: checked.
- No universal spectral-gap claim: checked.

