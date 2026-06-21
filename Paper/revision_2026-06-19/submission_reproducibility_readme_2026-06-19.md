# Submission reproducibility README — 2026-06-19

Scope: entry point for editors, reviewers, and future authors using the
submission materials for `Paper/revision_2026-06-19/draft.tex`.

This file is intentionally journal-neutral.  It describes how to navigate the
current reproducibility package before a final target journal, template, author
metadata, and possible DOI-backed raw-data archive are chosen.

## Current local status

The local submission gate currently reports:

```text
PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION
```

This means that the manuscript source, SIADS review-preparation source,
figures, bibliography, source-traced derived artifacts, claim audit,
availability-path audit, LaTeX compile/log scans, and minimal raw-data manifest
pass local checks.  It does **not** mean
that author declarations, professional similarity checking, journal-specific
formatting, or a DOI-backed raw-data archive have been completed.

## Fast verification path

From the repository root, run:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Expected current result:

- `latex_log`: `PASS`
- `siads_latex_log`: `PASS`
- `compiled_pdf_artifact_audit`: `PASS`
- `availability_path_audit`: `PASS`
- `manuscript_claim_audit`: `PASS`
- `reference_integrity_audit`: `PASS`
- `author_submission_fields_audit`: `AUTHOR_CONFIRMATION_PENDING`
- `raw_data_archive_manifest`: `PASS`
- `submission_bundle_manifest`: `PASS_WITH_LOCAL_RAW_DATA_LIMITATION`
- `submission_source_bundle`: `PASS`

The generated summary is:

- `Paper/revision_2026-06-19/submission_checks_summary.md`
- `Paper/revision_2026-06-19/submission_checks_summary.json`

## What each artifact is for

| Need | Start here | What it establishes |
|---|---|---|
| Read or compile the paper | `draft.tex`; `draft_siads_review.tex`; compiled PDFs from the LaTeX gates | Current manuscript text, SIADS review-preparation text, and figures. |
| Check compiled-PDF artifact metadata | `compiled_pdf_artifact_audit.md`; `compiled_pdf_artifact_audit.json` | The current generic and SIADS local PDFs have recorded page counts, byte sizes, and SHA-256 checksums. |
| Check the original working outline | `paper_draft_1.tex` | Preserved planning/framework source; not the active manuscript. |
| Check citation/reference integrity | `reference_integrity_audit.md`; `reference_integrity_audit.json` | All local cite keys and BibTeX entries match, identifiers are present, and external verification URLs are recorded. |
| Check numerical-claim support | `manuscript_claim_audit.md`; `manuscript_claim_audit.json` | Code-verifiable numerical/data claims currently pass 23/23 registered checks. |
| Check file/path availability | `availability_path_audit.md`; `availability_path_audit.json` | Manuscript-declared files and figure paths exist locally and have hashes where applicable. |
| Check submission-bundle completeness | `submission_bundle_manifest.md`; `submission_bundle_manifest.json` | Tracked release files are present, git-tracked, and categorized by role. |
| Check submission metadata consistency | `submission_metadata_consistency_audit.md`; `submission_metadata_consistency_audit.json` | Handoff documents quote the current PDF hashes, page counts, release-file count, and predicted source-bundle file count. |
| Build a source-only submission archive | `submission_source_bundle_report.md`; `scripts/build_submission_source_bundle.py` | A source-only `.tar.gz` can be generated under `tmp/`, with checksums and raw-data exclusions recorded. |
| Build a local journal-upload package | `journal_upload_file_index_2026-06-20.md`; `scripts/build_journal_upload_package.py` | A timestamped local package under `tmp/` collects the selected manuscript PDF, source archive, handoff documents, and checksum index for author/advisor upload review. |
| Build SIADS cover-letter template | `siads_cover_letter_template.tex`; `siads_cover_letter_template_build.md`; `scripts/build_siads_cover_letter_template.py` | A placeholder cover-letter PDF can be compiled locally for author editing; it is explicitly not a final submission letter. |
| Check compiled-PDF layout QA | `pdf_layout_qa_2026-06-19.md` | The generic compiled PDF has clean LaTeX logs and the newly added figure/table pages have been rendered and checked for obvious layout defects. |
| Check source-traced raw files | `raw_data_archive_manifest.md`; `raw_data_archive_manifest.json` | The compact raw-data subset contains 42 unique local files totaling 151,605,557 bytes. |
| Build a local raw-data upload archive | `scripts/build_raw_data_archive.py`; `raw_data_archive_build_report.md` | A timestamped `.tar.gz` preserving `raw_data/...` paths can be generated under `tmp/`; the latest local build archived 42/42 raw files with zero missing files. |
| Check author/journal-only blockers | `author_submission_fields_audit.md`; `author_submission_fields_audit.json` | Records the 9 current author/external items that must be resolved before formal submission. |
| Apply confirmed author fields | `author_submission_fields_template.json`; `scripts/apply_author_submission_fields.py` | Provides a dry-run-first workflow for replacing declaration/front-matter/data-availability placeholders after author confirmation. |
| Check current-scaling validation | `experiments/flux_validation/production_manifest.md`; `experiments/flux_validation/validation_report.md`; `experiments/flux_validation/production_dt5e-4/current_windows_window_statistics.csv` | Production flux/current scaling artifacts, validation summaries, and finite-window diagnostics. |
| Check gamma-robustness production path | `gamma_robustness_smoke_report.md`; `scripts/run_gamma_robustness_smoke.py` | Smoke-only verification that gamma-specific temporary sources can be generated from the frozen canonical source, compiled, and run without changing the production source hash. This is not manuscript evidence. |
| Check LTE residual mesh diagnostic | `report_assets/compare_residual_mesh.pdf`; `report_assets/compare_residual_mesh_metrics.md`; `draft.tex` | Structural LTE residual visualization including the `n=15` diagnostic requested in the manuscript revision pass, with source-traced mesh-slice RMS metrics. |
| Check author/journal blockers | `author_submission_action_packet_2026-06-19.md`; `target_journal_shortlist_2026-06-19.md`; `submission_readiness_checklist_2026-06-19.md` | Remaining human decisions before formal submission. |
| Prepare SIADS-first submission | `draft_siads_review.tex`; `siads_first_submission_packet_2026-06-20.md` | Line-numbered SIADS review-preparation source plus cover-letter draft, supplementary-material index, keywords/MSC candidates, and SIADS conversion checklist for the recommended first target. |
| See the full revision history | `progress_report.md`; `integrity_audit_2026-06-19.md` | What was changed, verified, limited, or left for authors. |

## Raw-data archive convention

If the final journal requires raw data beyond the tracked GitHub artifacts,
archive the files listed in `raw_data_archive_manifest.json` with their relative
paths preserved under:

```text
raw_data/
```

The Markdown manifest prints abbreviated SHA-256 digests for readability; the
JSON manifest records the full SHA-256 value for each raw file.  After creating
a Zenodo, OSF, institutional, or journal supplement archive, update the paper's
data/code availability statement with the archive DOI or accession link and
rerun the local gate.

The helper

```sh
python3 Paper/revision_2026-06-19/scripts/build_raw_data_archive.py
```

creates a timestamped local `.tar.gz` under `tmp/raw_data_archive/runs/` and
writes `raw_data_archive_build_report.md`.  This prepares the upload artifact;
it does not itself upload the archive or create a DOI.

## What is deliberately not claimed

The current package does not yet provide:

1. final author affiliations, ORCID IDs, corresponding-author metadata, funding
   statement, competing-interest statement, or final CRediT allocation;
2. professional plagiarism/self-plagiarism screening such as iThenticate or
   Turnitin;
3. a journal-specific LaTeX class conversion;
4. an uploaded DOI-backed raw-data archive;
5. full neural-network retraining reproducibility beyond saved-model inference
   and archived notebook/model outputs.

These items are tracked in
`author_submission_action_packet_2026-06-19.md` and
`submission_readiness_checklist_2026-06-19.md`.

## Suggested final handoff order

For a reviewer or coauthor who has not seen the project before:

1. read `draft.tex` or the compiled PDF;
2. read `submission_reproducibility_readme_2026-06-19.md` (this file);
3. run `scripts/run_submission_checks.py --compile-latex`;
4. inspect `submission_checks_summary.md`;
5. inspect `pdf_layout_qa_2026-06-19.md` for compiled-PDF layout status;
6. inspect `reference_integrity_audit.md` for citation/reference integrity;
7. inspect `submission_source_bundle_report.md` for source-archive packaging status;
8. inspect `manuscript_claim_audit.md` for claim-by-claim support;
9. inspect `submission_bundle_manifest.md` and
   `raw_data_archive_manifest.md` before deciding whether a raw-data supplement
   is needed;
10. complete the author-facing items in
   `author_submission_action_packet_2026-06-19.md`.
