# NLS final draft text + figures package

Prepared from `/Users/jayleenjiang/Documents/NLS/Paper/revision` on 2026-07-06.

This is the minimal package needed to compile the current main manuscript draft. It contains:

- `draft.tex`
- `references.bib`
- the 13 figure files directly referenced by `draft.tex`
- `compile_latex.sh`

It intentionally does not include raw simulation data, audit reports, auxiliary LaTeX files, or large experiment folders.

To compile:

```bash
./compile_latex.sh
```

or manually:

```bash
latexmk -pdf -interaction=nonstopmode draft.tex
```

If you compile manually with `pdflatex` and see warnings such as
`Reference ... undefined`, run the full LaTeX/BibTeX/LaTeX/LaTeX sequence or use
`latexmk`; cross-references are resolved on later passes.
