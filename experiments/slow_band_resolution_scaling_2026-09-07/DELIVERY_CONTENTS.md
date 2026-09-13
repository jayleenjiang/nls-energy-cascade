# Delivery contents

- `PROTOCOL.md`: pre-output scientific and interpretation rules.
- `run_scaling.py`: frozen Part-1 analysis implementation.
- `run_part1.sh`, `run_part1_when_ac.sh`: exact launch wrappers.
- `PART1_GATE.json`: binding Part-2 decision.
- `part1_analysis/synthetic_scaling_raw.csv`: all 168,000 replicate rows.
- `part1_analysis/synthetic_scaling_summary.csv`: all 336 cell summaries.
- `part1_analysis/delta020_summary.csv`: the 48 target-spacing cells.
- `part1_analysis/observed_resolution_passes.csv`: all 22 observed passing
  cells elsewhere on the tested grid.
- `part1_analysis/resolution_curve.csv`: 48 resolution-curve rows.
- `part1_analysis/delta020_extrapolation.csv`: eight descriptive scaling and
  saturation rows.
- `part1_analysis/task_seeds.csv`: all 84 deterministic task seeds.
- `part1_analysis/run_manifest.json`: command, hashes, grid, runtime, and gate.
- `audit_results.py`, `audit.log`: independent completeness and integrity
  audit.
- `provenance/integrity_audit.json`, `provenance/artifact_sha256.csv`: audit
  result and output hashes.
- `make_report.py`, `figures/`: report-generation source and all graphs in PDF
  and PNG formats.
- `report/slow_band_resolution_scaling_report.tex`: editable LaTeX report.
- `report/slow_band_resolution_scaling_report.pdf`: visually checked report.
- `FINAL_VERDICT.md`, `VALIDATION_REPORT.md`, `PROVENANCE.md`: claim boundary,
  audit summary, and reproducibility record.

The external delivery ZIP also contains the two hashed input NPZ correlation
arrays under `inputs/`; they are not duplicated inside this experiment folder.

