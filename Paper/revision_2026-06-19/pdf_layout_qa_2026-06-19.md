# PDF layout QA report — updated 2026-06-20

Scope: layout and rendering sanity check for the compiled PDFs generated from
`Paper/revision_2026-06-19/draft.tex` and
`Paper/revision_2026-06-19/draft_siads_review.tex`.  The visual-rendering
notes below are a targeted post-edit QA record after adding the LTE residual
decomposition table and the residual-mesh slice metrics; a final page-by-page
proof review is still required after author and journal edits.

This is a local production-quality check, not a target-journal proof review.
It should be rerun after any author metadata insertion, target-journal template
conversion, figure replacement, or final copyediting.

## Source and command

The PDF was generated from the repository root with:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

The latest runner reported
`PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`; the
LaTeX-specific gate reported `PASS`.

Current compiled PDFs after the LTE residual-mesh metrics update:

| PDF | Local build path | Pages | Size bytes | SHA-256 |
|---|---|---:|---:|---|
| Generic revised manuscript | `tmp/paper_build/revision/draft.pdf` | 27 | 1,572,459 | `88dcde6e7276cb903733e308d60ee4d3d8a8448893f24c5b7bb05c2c79800c29` |
| SIADS review-preparation source | `tmp/paper_build/siads_review/draft_siads_review.pdf` | 27 | 1,587,312 | `9dc6635c46dcddd4f8d45a584629d9c3878aba016ca4675c409654b601327f91` |

Both PDFs use A4 media boxes and were produced by pdfTeX.

## Machine checks

### LaTeX log scan

The generated log `tmp/paper_build/revision/draft.log` was scanned for:

- `LaTeX Warning`
- package warnings
- `Overfull`
- `Underfull`
- undefined references/citations

Result: no matches.

### TeX structural count

Counts from `draft.tex` after the LTE residual-decomposition update:

| Item | Count |
|---|---:|
| `figure` environments | 11 |
| `table` environments | 4 |
| `equation` environments | 18 |
| `align` environments | 1 |
| `includegraphics` calls | 12 |
| `\cite...` commands | 28 |

The 11 `figure` environments and 4 `table` environments all have captions in
the source.

## Visual rendering check

The current post-edit check focuses on the newly added LTE residual
even/odd-decomposition table.  The affected pages were located by PDF text
extraction, rendered from the compiled PDFs to PNG at 144 dpi with
`pdfplumber`, and inspected individually:

- generic manuscript, page 10: new table, caption, surrounding LTE text, and
  transition into the thermal-conductivity subsection;
- SIADS review-preparation PDF, page 10: same table with line numbers and the
  preceding LTE table on the same page.

Result:

- no page rotation problems;
- no visible black boxes after rendering with a white background;
- no visibly clipped newly added table or caption;
- no obvious overlapping text;
- the residual-decomposition table is readable in both the generic and SIADS
  review-preparation formats;
- the generic PDF page transition into `Action-current conductivity` is clean, and
  the SIADS page remains legible with line numbers.

The earlier 2026-06-19 full-page visual pass remains the baseline for the
unchanged pages; after any target-journal template conversion, rerun a full
page-by-page proof inspection.

## Remaining layout limitations

1. The PDF is not yet converted to a target-journal class or publisher template.
2. The PDF is not tagged for accessibility.
3. This check does not replace a final human proofread after journal-template
   conversion.
4. If the target journal imposes figure-resolution or color-space constraints,
   those constraints should be checked after selecting the journal.
