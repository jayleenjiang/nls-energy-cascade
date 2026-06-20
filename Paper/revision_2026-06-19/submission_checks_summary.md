# Submission checks summary

Generated: `2026-06-20T19:07:59.404206+00:00`

Overall status: **PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION**

## Gate summary

| Gate | Status | Key numbers |
|---|---|---|
| latex_log | PASS | issues=0 |
| availability_path_audit | PASS | total_paths=37, missing_paths=0, untracked_files=0 |
| manuscript_claim_audit | PASS | total_claims=18, verified=18, failed=0 |
| reference_integrity_audit | PASS | bib_entries=8, citation_commands=28, citation_uses=29, dangling_citation_count=0, missing_identifier_count=0, missing_required_field_count=0, missing_verification_source_count=0, orphan_reference_count=0, unique_cited_keys=8 |
| author_submission_fields_audit | AUTHOR_CONFIRMATION_PENDING | pass_count=0, pending_count=9, total_checks=9 |
| raw_data_archive_manifest | PASS | missing_file_count=0, referenced_file_count=40, referenced_total_bytes=138875181, status=PASS |
| submission_bundle_manifest | PASS_WITH_LOCAL_RAW_DATA_LIMITATION | local_raw_dependency_count=44, missing_release_files=0, release_directory_count=1, release_file_count=90, status=PASS_WITH_LOCAL_RAW_DATA_LIMITATION, untracked_local_raw_dependencies=44, untracked_release_files=0 |
| submission_source_bundle | PASS | directory_tracked_file_count=174, excluded_self_referential_count=6, included_file_count=259, included_total_bytes=4460866, local_raw_dependency_count=44, manifest_release_directories=1, manifest_release_file_records=90, missing_file_count=0 |

## Command results

| Command | Return code | Duration (s) |
|---|---:|---:|
| `latex_compile` | 0 | 0.061 |
| `availability_path_audit` | 0 | 0.178 |
| `manuscript_claim_audit` | 0 | 0.043 |
| `reference_integrity_audit` | 0 | 0.031 |
| `author_submission_fields_audit` | 0 | 0.029 |
| `submission_bundle_manifest_initial` | 0 | 0.786 |
| `raw_data_archive_manifest` | 0 | 0.186 |
| `submission_bundle_manifest_final` | 0 | 0.655 |
| `submission_source_bundle` | 0 | 0.325 |

## LaTeX log check

- Status: PASS

## Scope limitations

- This runner checks local reproducibility gates only.
- It does not confirm author/funding/competing-interest declarations.
- It does not run professional plagiarism/self-plagiarism screening.
- It does not upload raw data or create a DOI-backed archive.
