# Data-package contents

The complete local archive contains:

- this full experiment directory, including the four new 4,000,256-row raw
  block CSVs, burn-in/profile/summary outputs, failed-attempt and recovery
  logs, frozen protocol and amendments, analysis source, all raw-count and
  bootstrap tables, figures, TeX/PDF report, verdicts, and integrity audit;
- the three immutable `dt=5e-4` baseline block files used for weak, moderate,
  and equilibrium analysis;
- the immutable `dt=2.5e-4` equilibrium baseline block file;
- the exact production source and binary;
- the two frozen imported analysis modules and the completed `(7,5)` reference
  extrapolation table required by the analysis hash gates.

The archive preserves the repository-relative directory layout, so the
analysis command in `REPRODUCIBILITY_AUDIT.md` works after extraction at the
repository root.  LaTeX temporary files, Python bytecode, PID files, and run
locks are omitted; scientific logs and all immutable data are retained.
