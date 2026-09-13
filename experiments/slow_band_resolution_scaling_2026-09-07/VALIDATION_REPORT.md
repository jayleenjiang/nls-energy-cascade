# Validation report

Overall audit: **PASS**.

## Frozen-output checks

- **manifest complete**: PASS; status=complete.
- **summary cell count and grid**: PASS; rows=336, unique_cells=336, expected=336.
- **raw row count**: PASS; rows=168000, expected=168000.
- **raw replicate coverage**: PASS; cells=336, replicates_per_cell=500.
- **summary reconstructed from raw rows**: PASS; maximum_absolute_reconstruction_difference=0.000e+00.
- **finite summary statistics**: PASS; rows=336.
- **resolution curve row count**: PASS; rows=48, expected=48.
- **delta=0.20 extrapolation row count**: PASS; rows=8, expected=8.
- **frozen Part-2 gate reproduced**: PASS; delta020_passing_cells=0; part2_allowed=False.
- **analysis source hash**: PASS; sha256=b2ad9a90fc9ecfd2799c9b778dff3d7194c28faca224da70f20f5db022dce1d4.
- **previous pipeline hash**: PASS; sha256=ada94924e17f762519a973b28712a3af57ab6f9d1987183a2459c830587d94cb.
- **input NPZ hashes**: PASS; {"driven_dt2p5e-4": "9a878a637bd7475ab0f05ea0b7991276f4de1fd53f152c274530b7f99a984a47", "equilibrium_dt2p5e-4": "f7f97e860d75414b3999d656cf585dea9d402fafe7cc1a99a22321e9963c3d45"}.
- **deterministic task-seed coverage**: PASS; unique_case_N_separation_seeds=84, expected=84.

## Claim boundary

Part 2 was not run because the predeclared delta=0.20 gate failed. Power-law stream-count estimates are reported only as descriptive extrapolations and are not observed resolution results.
