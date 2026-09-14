# Reproducibility package contents

This directory contains the complete Section 4 analysis frozen on
2026-09-14.

- `PROTOCOL.md`: preregistered split, gate, and stop rules.
- `FINAL_VERDICT.md`, `VALIDATION_REPORT.md`: claim boundary and raw summaries.
- `paper/04-short-chain-rewritten.tex`: manuscript-ready rewritten section.
- `report/`: standalone TeX/PDF report and compile smoke wrapper.
- `figures/` and `analysis/`: publication figures and their plotted CSVs.
- `stationary_density/`: normalized density implementation, all candidate
  validation metrics/checkpoints, and the frozen negative selection result.
- `modes/`: frozen EDMD selections, one-shot blind outputs, correlations,
  matrices, and modal-weight bootstraps.
- `raw/direct_state_f64/`: all 256 float64 direct-state trajectories used by
  the analysis.
- `source/` and `bin/`: replay source and binary.
- `provenance/`: input/output hashes, exact commands, timestamps, and audit
  logs.

The large raw trajectory files are binary arrays of five float64 columns in
the order `(I1,I2,I3,theta1,theta3)`, with 100,001 rows per stream. Case
metadata and stream seeds are stored beside the arrays. `DIRECT_STATE_SPLITS.csv`
and `STREAM_SPLITS.csv` give the trajectory-level assignments.

The density result is intentionally negative. Do not train or select a NESS
density from the blind files in this package without declaring a new protocol.

