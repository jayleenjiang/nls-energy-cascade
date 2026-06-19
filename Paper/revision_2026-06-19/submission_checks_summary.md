# Submission checks summary

Generated: `2026-06-19T18:08:16.444524+00:00`

Overall status: **PASS_WITH_LOCAL_RAW_DATA_LIMITATION**

## Gate summary

| Gate | Status | Key numbers |
|---|---|---|
| latex_log | PASS | issues=0 |
| availability_path_audit | PASS | total_paths=31, missing_paths=0, untracked_files=0 |
| manuscript_claim_audit | PASS | total_claims=14, verified=14, failed=0 |
| raw_data_archive_manifest | PASS | missing_file_count=0, referenced_file_count=40, referenced_total_bytes=138875181, status=PASS |
| submission_bundle_manifest | PASS_WITH_LOCAL_RAW_DATA_LIMITATION | local_raw_dependency_count=44, missing_release_files=0, release_directory_count=1, release_file_count=53, status=PASS_WITH_LOCAL_RAW_DATA_LIMITATION, untracked_local_raw_dependencies=44, untracked_release_files=0 |

## Command results

| Command | Return code | Duration (s) |
|---|---:|---:|
| `latex_compile` | 0 | 0.063 |
| `availability_path_audit` | 0 | 0.168 |
| `manuscript_claim_audit` | 0 | 0.034 |
| `submission_bundle_manifest_initial` | 0 | 0.583 |
| `raw_data_archive_manifest` | 0 | 0.161 |
| `submission_bundle_manifest_final` | 0 | 0.517 |

## LaTeX log check

- Status: PASS

## Scope limitations

- This runner checks local reproducibility gates only.
- It does not confirm author/funding/competing-interest declarations.
- It does not run professional plagiarism/self-plagiarism screening.
- It does not upload raw data or create a DOI-backed archive.
