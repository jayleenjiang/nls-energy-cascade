# Submission reproducibility README — 2026-06-19

Scope: entry point for editors, reviewers, and future authors using the
submission materials for `Paper/revision_2026-06-19/draft.tex`.

This file is intentionally journal-neutral.  It describes how to navigate the
current reproducibility package before a final target journal, template, author
metadata, and possible DOI-backed raw-data archive are chosen.

## Current local status

The local submission gate currently reports:

```text
PASS_WITH_LOCAL_RAW_DATA_LIMITATION
```

This means that the manuscript source, figures, bibliography, source-traced
derived artifacts, claim audit, availability-path audit, LaTeX compile/log
scan, and minimal raw-data manifest pass local checks.  It does **not** mean
that author declarations, professional similarity checking, journal-specific
formatting, or a DOI-backed raw-data archive have been completed.

## Fast verification path

From the repository root, run:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Expected current result:

- `latex_log`: `PASS`
- `availability_path_audit`: `PASS`
- `manuscript_claim_audit`: `PASS`
- `reference_integrity_audit`: `PASS`
- `raw_data_archive_manifest`: `PASS`
- `submission_bundle_manifest`: `PASS_WITH_LOCAL_RAW_DATA_LIMITATION`
- `submission_source_bundle`: `PASS`

The generated summary is:

- `Paper/revision_2026-06-19/submission_checks_summary.md`
- `Paper/revision_2026-06-19/submission_checks_summary.json`

## What each artifact is for

| Need | Start here | What it establishes |
|---|---|---|
| Read or compile the paper | `draft.tex`; compiled PDF from the LaTeX gate | Current manuscript text and figures. |
| Check the original working outline | `paper_draft_1.tex` | Preserved planning/framework source; not the active manuscript. |
| Check citation/reference integrity | `reference_integrity_audit.md`; `reference_integrity_audit.json` | All local cite keys and BibTeX entries match, identifiers are present, and external verification URLs are recorded. |
| Check numerical-claim support | `manuscript_claim_audit.md`; `manuscript_claim_audit.json` | Code-verifiable numerical/data claims currently pass 18/18 registered checks. |
| Check file/path availability | `availability_path_audit.md`; `availability_path_audit.json` | Manuscript-declared files and figure paths exist locally and have hashes where applicable. |
| Check submission-bundle completeness | `submission_bundle_manifest.md`; `submission_bundle_manifest.json` | Tracked release files are present, git-tracked, and categorized by role. |
| Build a source-only submission archive | `submission_source_bundle_report.md`; `scripts/build_submission_source_bundle.py` | A source-only `.tar.gz` can be generated under `tmp/`, with checksums and raw-data exclusions recorded. |
| Check compiled-PDF layout QA | `pdf_layout_qa_2026-06-19.md` | The generic compiled PDF has clean LaTeX logs and the newly added figure/table pages have been rendered and checked for obvious layout defects. |
| Check source-traced raw files | `raw_data_archive_manifest.md`; `raw_data_archive_manifest.json` | The compact raw-data subset contains 40 unique local files totaling 138,875,181 bytes. |
| Build a local raw-data upload archive | `scripts/build_raw_data_archive.py`; `raw_data_archive_build_report.md` | A timestamped `.tar.gz` preserving `raw_data/...` paths can be generated under `tmp/`; the latest local build archived 40/40 raw files with zero missing files. |
| Check current-scaling validation | `experiments/flux_validation/production_manifest.md`; `experiments/flux_validation/validation_report.md`; `experiments/flux_validation/production_dt5e-4/current_windows_window_statistics.csv` | Production flux/current scaling artifacts, validation summaries, and finite-window diagnostics. |
| Check LTE residual mesh diagnostic | `report_assets/compare_residual_mesh.pdf`; `draft.tex` | Structural LTE residual visualization including the `n=15` diagnostic requested in the manuscript revision pass. |
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
