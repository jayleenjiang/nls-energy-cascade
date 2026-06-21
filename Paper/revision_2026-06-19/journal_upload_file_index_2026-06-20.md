# Journal upload file index — 2026-06-20

Scope: upload-facing index for the current near-submission revision of
`Paper/revision_2026-06-19/draft.tex` and the SIADS review-preparation source
`Paper/revision_2026-06-19/draft_siads_review.tex`.

This file does **not** mean the paper has been submitted.  It lists the local
files that are ready to hand to the authors/advisor for final review or upload
after author declarations, target-journal choices, and similarity screening are
completed.

## Current local gate

Latest command:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Latest status: `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`.

Key numbers:

| Gate item | Current value |
|---|---:|
| LaTeX/log checks | `draft.tex` PASS; `draft_siads_review.tex` PASS |
| Compiled PDF artifact audit | 2/2 PDFs verified |
| Availability paths | 38/38 |
| Registered numerical claims | 20/20 |
| Cited BibTeX entries | 8 |
| Release-bundle files | 121 |
| Source-bundle included files | 316 |
| Minimal raw-data files | 40 |
| Missing raw-data files | 0 |
| Author/submission-field audit | 9 pending author/external items |

## Manuscript PDFs

| Use | Local file | Pages | Size | SHA-256 | Status |
|---|---|---:|---:|---|---|
| Generic revised manuscript PDF | `tmp/paper_build/revision/draft.pdf` | 24 | 1,560,415 bytes | `bb20ce8e66dc99b83cf93026ea68f324ab89df24d9953e05c56ef670e8ff0ffe` | Suitable for author/advisor reading; recompile after final declarations. |
| SIADS review-preparation PDF | `tmp/paper_build/siads_review/draft_siads_review.pdf` | 24 | 1,575,302 bytes | `b37f974efc1389d8e63dbd5cb17c743ebafc150d954c7f362921f20ec6653b10` | Current SIADS-facing review PDF; use only after authors confirm SIADS and final metadata. |

The PDFs are local build artifacts under `tmp/` and are intentionally not
committed to Git.  Regenerate them after final author/journal edits.  The
machine-readable metadata for the current local PDFs is recorded in
`compiled_pdf_artifact_audit.json` and summarized in
`compiled_pdf_artifact_audit.md`.

## Source and reproducibility archives

| Use | Local file | Size | SHA-256 | Status |
|---|---|---:|---|---|
| Source-only submission bundle | See `submission_source_bundle_report.md` | See latest report | See latest report | Local upload-ready source archive; excludes large raw-data roots by design. The exact archive path/checksum changes whenever tracked handoff files change, so the report is the authoritative source. |
| Minimal raw-data archive | `tmp/raw_data_archive/runs/20260620T185610Z/NLS_numerical_study_raw_data_minimal.tar.gz` | 42,608,549 bytes | `1f8f8faa9bd9d73c804b51013549c63abc7af3a71febe05a31eb4df63ff4997f` | Local upload-ready raw-data archive; upload to Zenodo/OSF/journal storage only after authors choose this route. |

## Cover letter template

| Use | Local file | Status |
|---|---|---|
| SIADS cover-letter template PDF | `tmp/siads_cover_letter_template/siads_cover_letter_template.pdf` | Template only; contains placeholder fields and must be edited/approved before submission. |
| SIADS cover-letter template source | `Paper/revision_2026-06-19/siads_cover_letter_template.tex` | Tracked source template; compile status recorded in `siads_cover_letter_template_build.md`. |

The source bundle is documented by `submission_source_bundle_report.md`.  The
raw-data archive is documented by `raw_data_archive_build_report.md` and
preserves the `raw_data/...` paths listed in `raw_data_archive_manifest.md`.
For a local one-directory handoff, run
`scripts/build_journal_upload_package.py`; it creates a timestamped package
under `tmp/journal_upload_package/runs/` with the selected PDF, source archive,
handoff documents, and checksums.

## Files to upload by route

### Route A: SIADS review with repository-only data

Use only if the authors decide that a GitHub release is sufficient for the
initial review.

1. `tmp/paper_build/siads_review/draft_siads_review.pdf`
2. cover letter from `siads_first_submission_packet_2026-06-20.md`, after
   replacing bracketed author/funding/declaration fields
3. repository URL or GitHub release URL
4. optional source bundle:
   use the path and SHA-256 in `submission_source_bundle_report.md`

Do not cite a raw-data DOI in this route unless one has actually been created.

### Route B: SIADS review with DOI-backed raw-data supplement

Use if the authors or journal require archived raw data.

1. `tmp/paper_build/siads_review/draft_siads_review.pdf`
2. cover letter from `siads_first_submission_packet_2026-06-20.md`, after
   replacing bracketed author/funding/declaration fields
3. source bundle:
   use the path and SHA-256 in `submission_source_bundle_report.md`
4. upload raw-data archive to Zenodo/OSF/journal storage:
   `tmp/raw_data_archive/runs/20260620T185610Z/NLS_numerical_study_raw_data_minimal.tar.gz`
5. replace the manuscript data-availability statement with the resulting DOI
   or accession link, then rerun the full local gate

Do not upload the raw-data archive until the authors choose the archive route
and confirm any metadata required by the repository.

## Must be completed before formal upload

- Final author order, affiliations, corresponding author, email, and ORCID
  choices.
- Funding and competing-interest statements.
- Author-contribution approval wording.
- Professional similarity/self-plagiarism report.
- Target journal and article type.
- GitHub release tag and/or DOI-backed raw-data archive route.
- Final local gate after inserting all author/journal/release information.

The dry-run-first way to insert these fields is to fill
`author_submission_fields_template.json` as `author_submission_fields.json` and
then run `scripts/apply_author_submission_fields.py`.  Apply mode creates a
backup before editing both TeX sources.

## Regeneration commands

Run from the repository root:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
python3 Paper/revision_2026-06-19/scripts/build_submission_source_bundle.py
python3 Paper/revision_2026-06-19/scripts/build_raw_data_archive.py
python3 Paper/revision_2026-06-19/scripts/build_siads_cover_letter_template.py
python3 Paper/revision_2026-06-19/scripts/build_journal_upload_package.py --route siads-repository
```

The archive-building commands create timestamped artifacts under `tmp/`; update this
file if new archive paths or SHA-256 checksums are generated for final upload.
