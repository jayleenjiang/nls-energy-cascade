# Analysis-only roundoff-gate erratum

The frozen production simulation completed all eight batches and wrote all
`1,000,064` rows.  The first analysis invocation then stopped at its integrity
gate because recomputing
`Sigma_m=-Q_left/10-Q_right/2` from 17-digit decimal CSV fields differed from
the saved field by `2.842170943040401e-14`, while the script used a fixed
absolute threshold of `2e-14`.

This was an analysis implementation error: the threshold was smaller than the
roundoff accumulated when the three printed binary64 values were parsed and
combined.  The simulation, source, seed, parameters, and raw files were not
changed or rerun.  The complete failed output is preserved locally under
`analysis_failed_v1/` and in the tracked archive `analysis_failed_v1.zip`.

The repaired gate is scale-aware and fixed independently of the FT outcome:
`64 * machine_epsilon * max(1, magnitude of the compared fields)`.  The same
rule is now recorded for the balance identity and endpoint continuity checks.
All scientific support thresholds, histogram construction, bootstrap seeds,
and fit rules remain unchanged.  Only `analyze_feasibility.py` is rerun on the
original raw-file hashes.
