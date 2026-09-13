# NLS draft compile package

This package is a cleaned, self-contained bundle for compiling the main
manuscript `draft.tex`.

## Main files

- `draft.tex` — main manuscript source.
- `references.bib` — bibliography database used by `\bibliography{references}`.
- Figure files are stored at the same relative paths used in `draft.tex`.

## Compile

From this directory, run either:

```sh
latexmk -pdf -interaction=nonstopmode draft.tex
```

or:

```sh
pdflatex -interaction=nonstopmode draft.tex
bibtex draft
pdflatex -interaction=nonstopmode draft.tex
pdflatex -interaction=nonstopmode draft.tex
```

The package was prepared from `Paper/revision/` on 2026-07-06.

## Included material scope

Included:

- all images directly referenced by `draft.tex`;
- selected small CSV/JSON/MD support artifacts for current scaling, robustness,
  LTE residual metrics, canonical finite-time current distributions, validation
  checks, and audits;
- the latest local `submission_checks_summary.md`.

Not included:

- large raw-data roots such as `Energy Cascade/`, `KDE/`, and `lte/`;
- local binaries/build products not needed for compiling the manuscript;
- author-only final submission fields or external similarity-check reports.

The raw-data manifest in `raw_data_archive_manifest.md/json` records the compact
raw subset that would need a separate archive if a journal requires raw data.
