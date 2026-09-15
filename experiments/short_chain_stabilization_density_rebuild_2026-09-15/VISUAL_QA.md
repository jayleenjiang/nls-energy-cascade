# Visual QA

- `stabilization_density_rebuild_report.pdf`: five US-letter pages inspected
  after rendering with Poppler.  Tables and captions are readable; figures are
  ordered by section; no text, axes, legends, or color bars are clipped.
- `section4_2_smoke.pdf`: the standalone Section 4.2 fragment compiles and all
  three figures resolve.  Its last page contains the current figure by itself
  because the smoke wrapper has no surrounding manuscript text; this is normal
  float placement and not part of the fragment source layout contract.
- Both LaTeX logs contain no overfull boxes, underfull boxes, unresolved
  references, or undefined-control-sequence warnings after the final build.
- Figure PDFs use vector text and lines.  PNG copies are included for quick
  inspection only.
