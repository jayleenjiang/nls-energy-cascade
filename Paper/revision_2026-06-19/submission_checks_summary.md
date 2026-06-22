# Submission checks summary

Generated: `2026-06-22T04:47:27.739336+00:00`

Overall status: **PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION**

## Gate summary

| Gate | Status | Key numbers |
|---|---|---|
| latex_log | PASS | issues=0 |
| siads_latex_log | PASS | issues=0 |
| compiled_pdf_artifact_audit | PASS | pdf_count=2, failed=0 |
| availability_path_audit | PASS | total_paths=40, missing_paths=0, untracked_files=0 |
| manuscript_claim_audit | PASS | total_claims=23, verified=23, failed=0 |
| reference_integrity_audit | PASS | bib_entries=11, citation_commands=34, citation_uses=38, dangling_citation_count=0, missing_identifier_count=0, missing_required_field_count=0, missing_verification_source_count=0, orphan_reference_count=0, unique_cited_keys=11 |
| author_submission_fields_audit | AUTHOR_CONFIRMATION_PENDING | pass_count=0, pending_count=9, total_checks=9 |
| raw_data_archive_manifest | PASS | missing_file_count=0, referenced_file_count=42, referenced_total_bytes=151605557, status=PASS |
| submission_bundle_manifest | PASS_WITH_LOCAL_RAW_DATA_LIMITATION | local_raw_dependency_count=50, missing_release_files=0, release_directory_count=1, release_file_count=128, status=PASS_WITH_LOCAL_RAW_DATA_LIMITATION, untracked_local_raw_dependencies=50, untracked_release_files=0 |
| submission_metadata_consistency_audit | PASS | checked_documents=9, failed_checks=0, predicted_source_bundle_included_count=324, release_file_count=128, total_checks=21 |
| submission_source_bundle | PASS | directory_tracked_file_count=201, excluded_self_referential_count=6, included_file_count=324, included_total_bytes=6232585, local_raw_dependency_count=50, manifest_release_directories=1, manifest_release_file_records=128, missing_file_count=0 |
| siads_cover_letter_template | PASS | template_not_final=True, size_bytes=118094, issues=0 |
| journal_upload_package | PASS | contains_raw_data_archive=False, package_file_count=14, package_total_bytes=5370833 |

## Command results

| Command | Return code | Duration (s) |
|---|---:|---:|
| `latex_compile` | 0 | 0.058 |
| `siads_latex_compile` | 0 | 0.06 |
| `compiled_pdf_artifact_audit` | 0 | 0.049 |
| `availability_path_audit` | 0 | 0.185 |
| `manuscript_claim_audit` | 0 | 0.058 |
| `reference_integrity_audit` | 0 | 0.027 |
| `author_submission_fields_audit` | 0 | 0.023 |
| `submission_bundle_manifest_initial` | 0 | 0.955 |
| `raw_data_archive_manifest` | 0 | 0.099 |
| `submission_bundle_manifest_final` | 0 | 0.88 |
| `submission_metadata_consistency_audit` | 0 | 0.032 |
| `submission_bundle_manifest_post_metadata` | 0 | 0.908 |
| `submission_source_bundle` | 0 | 0.474 |
| `siads_cover_letter_template` | 0 | 0.076 |
| `journal_upload_package` | 0 | 0.125 |

## LaTeX log checks

- Generic manuscript: PASS
- SIADS review source: PASS

## Scope limitations

- This runner checks local reproducibility gates only.
- It does not confirm author/funding/competing-interest declarations.
- It does not run professional plagiarism/self-plagiarism screening.
- It does not upload raw data or create a DOI-backed archive.
