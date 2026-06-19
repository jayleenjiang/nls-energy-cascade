# PDF layout QA report — 2026-06-19

Scope: layout and rendering sanity check for the compiled PDF generated from
`Paper/revision_2026-06-19/draft.tex`.

This is a local production-quality check, not a target-journal proof review.
It should be rerun after any author metadata insertion, target-journal template
conversion, figure replacement, or final copyediting.

## Source and command

The PDF was generated from the repository root with:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

The runner reported `PASS_WITH_LOCAL_RAW_DATA_LIMITATION`; the LaTeX-specific
gate reported `PASS`.

Checked PDF:

- Local build path: `tmp/paper_build/revision/draft.pdf`
- SHA-256:
  `f8b2b0d8582730caf0648c3805374fbd48093628b73076f50355b17c0bad8e24`
- Size: 1,184,651 bytes
- Pages: 18
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

### PDF text-extraction sanity check

Text was extracted from all 18 pages with `pypdf`.

| Check | Count |
|---|---:|
| Extracted text characters | 48,612 |
| `??` markers | 0 |
| `TODO` | 0 |
| `TBD` | 0 |
| `FIXME` | 0 |
| `undefined` | 0 |
| `Missing` | 0 |
| `placeholder` / `PLACEHOLDER` | 0 |
| `References` heading | 1 |

The only capitalized `Reference` occurrence is the `References` heading on
page 18.

### TeX structural count

Counts from `draft.tex`:

| Item | Count |
|---|---:|
| `figure` environments | 9 |
| `table` environments | 2 |
| `equation` environments | 18 |
| `align` environments | 1 |
| `includegraphics` calls | 9 |
| `\cite...` commands | 28 |

The 9 `figure` environments and 2 `table` environments all have captions in
the source.

## Visual rendering check

All 18 pages were rendered to PNG at 110 dpi using the bundled Poppler
`pdftoppm` tool.  A contact sheet covering all pages was inspected, and the
most layout-sensitive pages were inspected individually:

- page 16: compact reproducibility table and beginning of data/code
  availability;
- page 17: availability list and declarations;
- page 18: references.

Result:

- no blank pages;
- no page rotation problems;
- no visible black boxes or missing glyph blocks;
- no visibly clipped figures or tables;
- no obvious overlapping text;
- no visibly truncated references;
- the compact reproducibility table on page 16 is dense but readable in the
  generic `article` format.

## Remaining layout limitations

1. The PDF is not yet converted to a target-journal class or publisher template.
2. The PDF is not tagged for accessibility.
3. This check does not replace a final human proofread after journal-template
   conversion.
4. If the target journal imposes figure-resolution or color-space constraints,
   those constraints should be checked after selecting the journal.
