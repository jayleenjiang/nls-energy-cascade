# Report-formatting amendment

Date: 2026-09-12.

After all frozen simulation and analysis outputs were complete, visual PDF
inspection found that the run-level provenance table began near the foot of
page 8 and continued on page 9 without a repeated header.  A `\clearpage` was
inserted immediately before that section so the complete table and provenance
text render together.  No simulation output, analysis value, statistical
gate, interpretation, figure, or frozen scientific parameter was changed.
