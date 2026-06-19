# Submission checks summary

Generated: `2026-06-19T18:13:40.114081+00:00`

Overall status: **PASS_WITH_LOCAL_RAW_DATA_LIMITATION**

## Gate summary

| Gate | Status | Key numbers |
|---|---|---|
| latex_log | PASS | issues=0 |
| availability_path_audit | PASS | total_paths=31, missing_paths=0, untracked_files=0 |
| manuscript_claim_audit | PASS | total_claims=14, verified=14, failed=0 |
| reference_integrity_audit | PASS | bib_entries=8, citation_commands=28, citation_uses=29, dangling_citation_count=0, missing_identifier_count=0, missing_required_field_count=0, missing_verification_source_count=0, orphan_reference_count=0, unique_cited_keys=8 |
| raw_data_archive_manifest | PASS | missing_file_count=0, referenced_file_count=40, referenced_total_bytes=138875181, status=PASS |
| submission_bundle_manifest | PASS_WITH_LOCAL_RAW_DATA_LIMITATION | local_raw_dependency_count=44, missing_release_files=0, release_directory_count=1, release_file_count=56, status=PASS_WITH_LOCAL_RAW_DATA_LIMITATION, untracked_local_raw_dependencies=44, untracked_release_files=0 |

## Command results

| Command | Return code | Duration (s) |
|---|---:|---:|
| `latex_compile` | 0 | 0.06 |
| `availability_path_audit` | 0 | 0.171 |
| `manuscript_claim_audit` | 0 | 0.038 |
| `reference_integrity_audit` | 0 | 0.028 |
| `submission_bundle_manifest_initial` | 0 | 0.656 |
| `raw_data_archive_manifest` | 0 | 0.158 |
| `submission_bundle_manifest_final` | 0 | 0.525 |

## LaTeX log check

- Status: PASS

## Scope limitations

- This runner checks local reproducibility gates only.
- It does not confirm author/funding/competing-interest declarations.
- It does not run professional plagiarism/self-plagiarism screening.
- It does not upload raw data or create a DOI-backed archive.
