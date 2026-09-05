# Advisor-safe migration draft

This directory is a self-contained compile package based on the advisor-synced
file `draft (1).tex`. The source file `main (7).tex` was used only as a source
of candidate prose and organization. No original file was overwritten.

## Main files

- `draft_advisor_safe.tex`: migrated manuscript source.
- `draft_advisor_safe.pdf`: verified 15-page build.
- `Figures/`: every image required by the TeX source.
- `MIGRATION_LOG.md`: included, excluded, and still-pending material.
- `source_snapshots/`: unchanged snapshots of both input manuscripts and the
  MATLAB residual scripts used to generate the LTE mesh.

## Compile

From this directory, run:

```bash
latexmk -pdf -interaction=nonstopmode -halt-on-error draft_advisor_safe.tex
```

The verified build used TeX Live 2023 and completed without undefined
citations, undefined cross-references, overfull boxes, or missing figures.

