# Validation report: one-bath open-chain profile experiment

## Production integrity

- Frozen remote implementation commit:
  `f4b74a06e658478e78b130975b34af20da9d9398`.
- Open-chain source SHA-256:
  `dfee7ed0b416a30f2b27d919c524dfaf37e91ca44312ffbb9d68b7fc1cd36110`.
- Binary SHA-256:
  `09d801e71b7f3bd13cea56aba3d432b471dbf16a726a9866894148a91397575d`.
- Twelve frozen replicate runs completed: two temperatures, three chain
  lengths, and two seeds per logical condition.
- All 3,072 requested trajectories are valid; non-finite trajectories: 0;
  discarded trajectories: 0.
- Every run used two threads, `gamma=0.1`, maximum `dt=5e-4`, measurement time
  2000, and the frozen size-dependent burn-in.

The mechanical audit in `audit_outputs.py` checks all 108 raw CSVs.  Every file
has its exact predeclared header and row count.  Trajectory IDs are contiguous
from 0 through 255 in every replicate; every profile/checkpoint sample count is
256.  The only blank fields are the two bond-sine fields at `j=n`, where no
bond exists.  All 108 hashes match both `RAW_HASHES.csv` and the independent
analysis input manifest.  The audit result is `PASS` with an empty error list
in `analysis/VALIDATION_REPORT.json`.

The detailed raw outputs record the explicit boundary equations, including
`right:b_I=0`, `right:b_phi=0`, `right:sigma_I=0`,
`right:sigma_phi=0`, and `right_boundary_bath=none`.  The summary rows record
boundary ID 4 and label `OPEN_left_BC1_right_free`.

## Analysis checks

- Six merged logical profiles and 30 curated CSVs were produced.
- Separate Q1--Q4 profile averages, instantaneous burn-in/measurement
  snapshots, total action, endpoint actions, and paired Q4--Q3 differences are
  present for every condition.
- The frozen stationarity gate was applied without exclusions or threshold
  changes.  Five conditions pass and `T10_n100` fails only the free-end CI-zero
  subgate (`z=-2.077`).
- Scaling uses exactly `n=25,50,100`; each exponent fit is explicitly marked as
  having one residual degree of freedom.
- Collapse uses the fixed exponent 1/2 and fixed grid `x=0,0.01,...,1`.
- Penetration uses the predeclared 50% threshold without extrapolation.

## Report and package checks

The nine-page PDF was rendered to PNG and visually inspected page by page.
Plots, confidence bands, equations, tables, labels, and hashes are legible and
not clipped.  A formatting-only amendment moved the complete provenance table
to a fresh page; no scientific output or gate changed.  The ZIP archive passes
an integrity test with no compressed-data errors.
