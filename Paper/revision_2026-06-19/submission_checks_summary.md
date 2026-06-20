# Submission checks summary

Generated: `2026-06-20T19:53:43.996975+00:00`

Overall status: **PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION**

## Gate summary

| Gate | Status | Key numbers |
|---|---|---|
| latex_log | PASS | issues=0 |
| siads_latex_log | PASS | issues=0 |
| compiled_pdf_artifact_audit | PASS | pdf_count=2, failed=0 |
| availability_path_audit | PASS | total_paths=37, missing_paths=0, untracked_files=0 |
| manuscript_claim_audit | PASS | total_claims=19, verified=19, failed=0 |
| reference_integrity_audit | PASS | bib_entries=8, citation_commands=28, citation_uses=29, dangling_citation_count=0, missing_identifier_count=0, missing_required_field_count=0, missing_verification_source_count=0, orphan_reference_count=0, unique_cited_keys=8 |
| author_submission_fields_audit | AUTHOR_CONFIRMATION_PENDING | pass_count=0, pending_count=9, total_checks=9 |
| raw_data_archive_manifest | PASS | missing_file_count=0, referenced_file_count=40, referenced_total_bytes=138875181, status=PASS |
| submission_bundle_manifest | PASS_WITH_LOCAL_RAW_DATA_LIMITATION | local_raw_dependency_count=44, missing_release_files=0, release_directory_count=1, release_file_count=95, status=PASS_WITH_LOCAL_RAW_DATA_LIMITATION, untracked_local_raw_dependencies=44, untracked_release_files=0 |
| submission_source_bundle | PASS | directory_tracked_file_count=174, excluded_self_referential_count=6, included_file_count=264, included_total_bytes=4524344, local_raw_dependency_count=44, manifest_release_directories=1, manifest_release_file_records=95, missing_file_count=0 |

## Command results

| Command | Return code | Duration (s) |
|---|---:|---:|
| `latex_compile` | 0 | 0.079 |
| `siads_latex_compile` | 0 | 0.057 |
| `compiled_pdf_artifact_audit` | 0 | 0.036 |
| `availability_path_audit` | 0 | 0.162 |
| `manuscript_claim_audit` | 0 | 0.046 |
| `reference_integrity_audit` | 0 | 0.028 |
| `author_submission_fields_audit` | 0 | 0.023 |
| `submission_bundle_manifest_initial` | 0 | 1.086 |
| `raw_data_archive_manifest` | 0 | 0.184 |
| `submission_bundle_manifest_final` | 0 | 0.745 |
| `submission_source_bundle` | 0 | 0.328 |

## LaTeX log checks

- Generic manuscript: PASS
- SIADS review source: PASS

## Scope limitations

- This runner checks local reproducibility gates only.
- It does not confirm author/funding/competing-interest declarations.
- It does not run professional plagiarism/self-plagiarism screening.
- It does not upload raw data or create a DOI-backed archive.
