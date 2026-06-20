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

Latest status: `PASS_WITH_LOCAL_RAW_DATA_LIMITATION`.

Key numbers:

| Gate item | Current value |
|---|---:|
| Availability paths | 37/37 |
| Registered numerical claims | 18/18 |
| Cited BibTeX entries | 8 |
| Release-bundle files | 87 |
| Source-bundle included files | 256 |
| Minimal raw-data files | 40 |
| Missing raw-data files | 0 |

## Manuscript PDFs

| Use | Local file | Pages | Size | SHA-256 | Status |
|---|---|---:|---:|---|---|
| Generic revised manuscript PDF | `tmp/paper_build/revision/draft.pdf` | 21 | not committed; local build artifact | `be3ce88e36f3b2324eeb058e8d2b5e9a74f072e9555cf5dd82e3d7a9ec45a37c` | Suitable for author/advisor reading; recompile after final declarations. |
| SIADS review-preparation PDF | `tmp/paper_build/siads_review/draft_siads_review.pdf` | 21 | about 1.5 MB | `1534e68e0617d0f44bc3d4b856d1ee608029ba3cc84c409a64802cfc93a4950b` | Current SIADS-facing review PDF; use only after authors confirm SIADS and final metadata. |

The PDFs are local build artifacts under `tmp/` and are intentionally not
committed to Git.  Regenerate them after final author/journal edits.

## Source and reproducibility archives

| Use | Local file | Size | SHA-256 | Status |
|---|---|---:|---|---|
| Source-only submission bundle | See `submission_source_bundle_report.md` | See latest report | See latest report | Local upload-ready source archive; excludes large raw-data roots by design. The exact archive path/checksum changes whenever tracked handoff files change, so the report is the authoritative source. |
| Minimal raw-data archive | `tmp/raw_data_archive/runs/20260620T185610Z/NLS_numerical_study_raw_data_minimal.tar.gz` | 42,608,549 bytes | `1f8f8faa9bd9d73c804b51013549c63abc7af3a71febe05a31eb4df63ff4997f` | Local upload-ready raw-data archive; upload to Zenodo/OSF/journal storage only after authors choose this route. |

The source bundle is documented by `submission_source_bundle_report.md`.  The
raw-data archive is documented by `raw_data_archive_build_report.md` and
preserves the `raw_data/...` paths listed in `raw_data_archive_manifest.md`.

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

## Regeneration commands

Run from the repository root:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
python3 Paper/revision_2026-06-19/scripts/build_submission_source_bundle.py
python3 Paper/revision_2026-06-19/scripts/build_raw_data_archive.py
```

The last two commands create timestamped archives under `tmp/`; update this
file if new archive paths or SHA-256 checksums are generated for final upload.
