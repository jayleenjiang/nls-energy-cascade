#!/usr/bin/env bash
set -euo pipefail

cd "$(dirname "$0")"

if command -v latexmk >/dev/null 2>&1; then
  latexmk -pdf -interaction=nonstopmode draft.tex
else
  pdflatex -interaction=nonstopmode draft.tex
  bibtex draft
  pdflatex -interaction=nonstopmode draft.tex
  pdflatex -interaction=nonstopmode draft.tex
fi
