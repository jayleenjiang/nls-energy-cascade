# Curated entropy/action fluctuation results

This directory contains the GitHub-sized, audited products of the long
`T_left=10`, `T_right=2` run.  The numerical data have passed the production
audits, but the study remains **provisional with respect to any fluctuation-
theorem claim**: the accessible finite-time data do not verify the exact FT.
The four raw `*_blocks.csv` files remain local and are identified by
`raw_blocks_sha256.txt`; they are not duplicated here.

The accepted local analysis package is
`../production/final_v4/`.  `final_v1` is a preserved failed layout gate,
while `final_v2` and `final_v3` are successful intermediate finalizations.
No scientific failure occurred in those earlier directories.

Key files:

- `MENTOR_SUMMARY.md`: meeting-ready interpretation and scope limits;
- `PAPER_SECTION_CANDIDATE.tex`: a standalone candidate subsection, not yet
  merged into the mentor-synchronized manuscript;
- `entropy_ft_report.pdf`: full 19-page audited report;
- `production_means.csv`: compact mean-current/entropy table;
- `ft_summary.csv`, `heat_symmetry_summary.csv`, and
  `action_symmetry_summary.csv`: primary fixed-bin symmetry estimates;
- `adaptive_symmetry_summary.csv`: independent fit-range robustness estimates;
- `action_tail_time_scaling.csv`: descriptive fits of
  `log P = intercept - t I(A)`;
- `action_normal_tail_fit_metrics.csv`: full-sample and joint two-tail Gaussian
  benchmarks;
- `coupling_summary.csv` and `stationarity_summary.csv`: heat--action coupling
  and drift checks;
- `raw_audit.md` and `analysis_audit.md`: fail-fast audit verdicts.

Interpretation rule: the medium-entropy and heat-current slopes below their
simple asymptotic references are reported as **not verified in the sampled
window**.  The experiment does not reconstruct the NESS system-entropy
endpoint term needed for an exact finite-time total-entropy FT.
