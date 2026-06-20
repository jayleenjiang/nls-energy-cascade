# PDF layout QA report — updated 2026-06-20

Scope: layout and rendering sanity check for the compiled PDF generated from
`Paper/revision_2026-06-19/draft.tex`.  The visual-rendering notes below are
a baseline layout QA record; the current 2026-06-20 interpretation edit is
text-only and still requires a final page-by-page proof review after author
and journal edits.

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

Current compiled PDF after the latest text-only edit:

- Local build path: `tmp/paper_build/revision/draft.pdf`
- SHA-256:
  `8c729b647149ee4df690f1c76d566909ce861c1ff98a9f40764728522df3768a`
- Size: 1,546,086 bytes
- Pages: 21
- Page size: A4
- PDF version: 1.5
- Producer: `pdfTeX-1.40.25`

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

Counts from `draft.tex`:

| Item | Count |
|---|---:|
| `figure` environments | 11 |
| `table` environments | 3 |
| `equation` environments | 18 |
| `align` environments | 1 |
| `includegraphics` calls | 12 |
| `\cite...` commands | 28 |

The 11 `figure` environments and 3 `table` environments all have captions in
the source.

## Visual rendering check

The updated manuscript adds the LTE residual mesh figure, the timestep
sensitivity table, and the finite-window current diagnostics figure.  These new
layout-sensitive pages were rendered from the compiled PDF to PNG at 144 dpi
with ImageMagick and inspected individually:

- page 11: LTE residual mesh figure and timestep sensitivity table;
- page 13: finite-window current diagnostics figure and short-chain section
  transition.

Result:

- no page rotation problems;
- no visible black boxes after rendering with a white background;
- no visibly clipped newly added figures or tables;
- no obvious overlapping text;
- the newly added timestep sensitivity table is dense but readable in the
  generic `article` format.

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
