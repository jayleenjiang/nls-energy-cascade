# Package contents

The package contains the frozen protocol, post-hoc current-scope audit,
analysis and audit source, exact commands, all recomputed CSV/JSON tables,
publication figures, the detailed report (`.tex` and `.pdf`), and the
paper-ready Section 4.2 fragment with a standalone smoke-test PDF.

The multi-gigabyte held-out trajectories and trained model files are not
duplicated.  Their exact paths, byte counts, and SHA-256 values are recorded in
`analysis/input_hashes.csv`; the V6 split manifest and density manifest hashes
are included there as well.

`analysis/output_integrity.json` records row counts and finite-value checks.
`analysis/artifact_hashes.csv` hashes every packaged analysis, figure, source,
and report artifact before packaging.  The ZIP SHA-256 is written beside the
archive.
