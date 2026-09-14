# Analysis erratum 002: omitted blind modal weights

Date: 2026-09-14

The one-shot Section-4.3 blind evaluator completed all frozen residual,
autocorrelation, plateau, damped-oscillation, AIC, bootstrap, and timestep
tests, but the implementation omitted the four modal-weight diagnostics named
in the protocol.  `MODE_SELECTION.json`, its SHA-256, the test results, and the
original completed read marker are preserved unchanged.

Before reopening any test file, this erratum freezes one additional
deterministic pass.  It may compute only the modal contributions of the
already selected modes to the already declared observables

1. `I2`,
2. `I1+I3`,
3. `cos(theta3)`, and
4. `cos(theta1)+cos(theta3)`.

The selected dictionary, cutoff, lag, preprocessing, eigenvalue, and
eigenvector remain fixed.  The pass uses the same EDMD modal-decomposition
formula as `edmd_oscillatory_extension_2026-09-07/run_extended_edmd.py`:
the test covariance and observable cross-moments are centered, transformed to
the frozen EDMD eigenbasis, and the selected eigenmode contribution (including
its conjugate partner for a complex mode) is divided by the EDMD-reconstructed
observable variance.  Uncertainty is a 500-replicate whole-test-stream
bootstrap with seed 2026091405.

This additional read cannot select a mode, change a gate, rescue a failed
mode, or authorize a unique spectral-gap claim.  It exists solely to complete
the predeclared diagnostic table.  The extra read is recorded separately in
`modes/MODAL_WEIGHT_READ_MARKER.json`.

