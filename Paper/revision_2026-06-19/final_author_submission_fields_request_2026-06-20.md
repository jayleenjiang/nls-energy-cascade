# Final author/journal information request — 2026-06-20

Scope: concise confirmation form for turning the locally verified manuscript
`Paper/revision_2026-06-19/draft.tex` and synchronized
`draft_siads_review.tex` into a formally submittable journal manuscript.

The local numerical, reference, LaTeX, path, source-bundle, and raw-data
manifest checks are already scripted.  The remaining items below require
author, advisor, journal-system, or external-service confirmation.  Do not
replace the manuscript declarations or set any confirmation booleans to `true`
until the corresponding action has actually been completed.

## Current local gate

Run from the repository root:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Current expected status before author/external items are resolved:

```text
PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION
```

This status means the local scientific/reproducibility gate is passing, but it
does not certify author approval, journal formatting, professional similarity
screening, or DOI-backed raw-data upload.

## Information to confirm

Fill these items first in prose.  After all answers are confirmed, copy
`author_submission_fields_template.json` to `author_submission_fields.json` and
enter the corresponding values there.

### 1. Target journal

- Target journal:
- Article type:
- Final format route:
  - generic article source for initial advisor circulation;
  - SIADS line-numbered review source;
  - official journal class/template conversion;
  - other:
- If SIADS is chosen, should the first upload use the current line-numbered
  fallback source or wait for conversion to official SIAM macros?
- Does the journal require anonymized review? yes / no / unclear
- Does the journal require ORCID identifiers? yes / no / unclear
- Does the journal require a data/code DOI? yes / no / unclear

JSON fields this controls:

- `target_journal`
- `article_type`
- `target_journal_confirmed`
- `notes`

### 2. Author metadata and approval

- Final author order:
- Jayleen Jiang affiliation:
- Jayleen Jiang email:
- Jayleen Jiang ORCID, if any:
- Yao Li affiliation:
- Yao Li email:
- Yao Li ORCID, if any:
- Corresponding author:
- Corresponding-author email:
- Have all authors reviewed and approved the final submitted version? yes / no

JSON fields this controls:

- `author_latex`
- `author_contributions_tex`
- `author_approval_confirmed`
- `notes`

Current recommended contribution wording, if accurate:

```text
Jayleen Jiang performed the numerical experiments, assembled the computational
artifacts, and drafted the manuscript. Yao Li supervised the project and
contributed to the model formulation, theoretical framing, and interpretation.
Both authors reviewed and approved the final submitted version.
```

### 3. Declarations

- Funding:
  - no external funding;
  - or grant/funder names and award numbers:
- Competing interests:
  - none;
  - or disclosure text:
- Ethics statement:
  - confirm current numerical/theoretical-study statement;
  - or supply target-journal wording:
- AI-assisted preparation statement:
  - confirm current wording;
  - or supply target-journal wording:

JSON fields this controls:

- `funding_tex`
- `competing_interests_tex`
- `notes`

Suggested final wording only if factually accurate:

```text
The authors declare no competing interests.
```

```text
The authors received no external funding for this work.
```

### 4. Data/code release route

Choose exactly one route before final submission.

| Route | Use when | Required final fields |
|---|---|---|
| `github_release` | The journal accepts a GitHub release/tag plus tracked reproducibility artifacts. | repository release URL |
| `doi_archive` | The journal expects DOI-backed raw data or stronger archival availability. | repository release URL, raw-data DOI/accession, upload confirmation |

The local minimal raw-data archive manifest currently covers 40 unique
source-trace raw files totaling 138,875,181 bytes.  A local upload-ready archive
can be built with:

```sh
python3 Paper/revision_2026-06-19/scripts/build_raw_data_archive.py
```

The archive must still be uploaded to Zenodo, OSF, or another approved service
before a DOI/accession can be cited.

JSON fields this controls:

- `data_release_route`
- `data_availability_tex`
- `repository_release_url`
- `raw_data_doi_or_accession`
- `raw_data_upload_completed_if_applicable`

Do not set `raw_data_upload_completed_if_applicable` to `true` unless the
chosen route is `doi_archive` and the external upload is complete.

### 5. Professional similarity/self-plagiarism screening

- Service used, e.g. iThenticate/Turnitin/journal system:
- Date completed:
- Result reviewed by:
- Any required changes made? yes / no / not applicable

JSON field this controls:

- `similarity_screening_completed`

The local originality spot-check is only a pre-screen and does not replace this
professional check.

### 6. Final PDF review

After all declaration, author, data-release, and journal-format edits:

- Recompile both manuscript sources.
- Review the final PDF page by page.
- Confirm figures, tables, captions, equations, references, line numbers, and
  declarations render correctly.

JSON field this controls:

- `final_pdf_review_completed`

## Dry-run and apply workflow

After author/journal information is confirmed:

```sh
cp Paper/revision_2026-06-19/author_submission_fields_template.json \
   Paper/revision_2026-06-19/author_submission_fields.json
```

Edit `author_submission_fields.json`, then run the dry run:

```sh
python3 Paper/revision_2026-06-19/scripts/apply_author_submission_fields.py \
  --input Paper/revision_2026-06-19/author_submission_fields.json
```

Apply only after the dry run reports `ready_to_apply=true`:

```sh
python3 Paper/revision_2026-06-19/scripts/apply_author_submission_fields.py \
  --input Paper/revision_2026-06-19/author_submission_fields.json \
  --apply
```

The apply command backs up `draft.tex` and `draft_siads_review.tex` before
editing.

Then rerun the full local gate:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Formal submission should proceed only after the gate still passes, the
author/submission-field audit no longer reports placeholder declarations, and
the final compiled PDF has been visually checked.
