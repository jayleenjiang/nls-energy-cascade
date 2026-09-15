## Material Passport

- Artifact type: numerical-experiment validation report
- Status: VERIFIED
- Data access: raw local float32 trajectories plus derived CSV/JSON
- Frozen remote commit: `fddd8c5a1700675a702aadab5dac9262ab48f6c6`

# Validation report

- Verification status: **VERIFIED** (production completed; integrity and
  reproducibility audits passed).
- Frozen remote protocol/source commit:
  `fddd8c5a1700675a702aadab5dac9262ab48f6c6`.
- Four cases completed: two bath conditions by two timesteps, each with 64
  streams and 100001 snapshots per stream.
- Projection/floor/zero-radius/midpoint-failure/nonfinite counts: all zero.
- All stored values finite; binary shapes and stream IDs pass.
- Worst discard-first-quarter stationarity metric: `0.00498164`, below the
  frozen `0.05` gate.
- Smallest effective count at an accepted plateau's terminal lag: `53259`.
- Common plateau: **FAIL** in all four cases under the unchanged rule.
- Timestep consistency: **PASS** for every observable resolved at both steps.
- Early `I1-I3` oscillation: **CONFIRMED**, with frequency close to the
  historical 5.35 at both steps and both bath conditions.
- Claim boundary: observable-specific rates and oscillatory modes are resolved;
  a unique numerical spectral gap is not.

All raw stream hashes and detailed case audits are in
`provenance/raw_artifact_manifest.sha256` and
`provenance/integrity_audit.json`.

## Statistical fallacy scan

1. **P-value-only inference:** not used; rates, windows, and bootstrap CIs are
   reported.
2. **Effect-size omission:** not applicable; decay rates and their differences
   are the reported effects.
3. **Correlation as causation:** avoided; the experiment distinguishes
   numerical sensitivity from persistence under a controlled scheme, not a
   formal causal proof about every historical implementation detail.
4. **Multiple comparisons:** no multiplicity-adjusted discovery claim is made;
   the predeclared common-rate gate requires simultaneous CI intersection.
5. **Post-selection/tuned windows:** avoided; plateau and damped-fit windows
   were frozen before production.
6. **Optional stopping:** avoided; trajectory count and duration were frozen.
7. **Pseudoreplication:** avoided in uncertainty estimates; bootstrap units are
   the 64 independent trajectories, not individual snapshots.
8. **Ignored serial correlation:** addressed by trajectory-level bootstrap and
   integrated-autocorrelation effective counts.
9. **Non-rejection as equality:** avoided; CI overlap plus a 20% numerical
   tolerance is reported as timestep consistency, not exact equality.
10. **Unsupported extrapolation:** avoided; no late-time extrapolation beyond
    directly supported lags is used to quote a gap.
11. **Scope/generalization error:** avoided; conclusions are restricted to the
    two bath cases, two timesteps, and accessible lag window for `n=3`.
