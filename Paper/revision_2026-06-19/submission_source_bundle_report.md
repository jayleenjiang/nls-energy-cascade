# Submission source bundle report

Generated: `2026-06-20T13:00:15.943552+00:00`

Status: **PASS**

## Archive

- Archive path: `tmp/submission_source_bundle/runs/20260620T130015Z/NLS_numerical_study_source.tar.gz`
- Archive bytes: 2470739
- Archive SHA-256: `30f83aa958007b5b8c64cee2f22b6b9541ca081a1dd12187c6f0abd12e26bab7`
- Staging directory: `tmp/submission_source_bundle/runs/20260620T130015Z/NLS_numerical_study_source`

## Summary

| Metric | Count |
|---|---:|
| Manifest release-file records | 64 |
| Manifest release directories | 1 |
| Included regular files | 172 |
| Included bytes | 3559239 |
| Excluded volatile/self-referential files | 6 |
| Missing files | 0 |
| Directory-tracked files copied | 113 |
| Local raw dependency records excluded | 44 |

## Excluded self-referential files

- `Paper/revision_2026-06-19/submission_bundle_manifest.json`
- `Paper/revision_2026-06-19/submission_bundle_manifest.md`
- `Paper/revision_2026-06-19/submission_checks_summary.json`
- `Paper/revision_2026-06-19/submission_checks_summary.md`
- `Paper/revision_2026-06-19/submission_source_bundle_report.json`
- `Paper/revision_2026-06-19/submission_source_bundle_report.md`

## Missing files

- None.

## Notes

- The archive is source-only and intentionally excludes the large local raw-data roots.
- If a target journal requires raw data, use `raw_data_archive_manifest.md` to prepare a DOI-backed raw-data supplement.
- The archive is written under `tmp/` and is not intended to be committed to Git.
- The file `SOURCE_BUNDLE_CONTENTS.sha256` inside the archive records SHA-256 checksums for copied files.
