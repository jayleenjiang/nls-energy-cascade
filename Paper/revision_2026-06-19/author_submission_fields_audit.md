# Author/submission field audit

Generated: `2026-06-21T02:46:36.821759+00:00`

Status: **AUTHOR_CONFIRMATION_PENDING**

## Summary

| Metric | Count |
|---|---:|
| Checks | 9 |
| Pending author/external items | 9 |
| Passed items | 0 |

## Checks

| ID | Category | Status | Evidence | Required action |
|---|---|---|---|---|
| `author_approval_placeholder` | manuscript-declarations | PENDING_AUTHOR_OR_EXTERNAL_ACTION | Provisional approval sentence present in draft.tex, draft_siads_review.tex. | After both authors approve, replace with final target-journal wording. |
| `competing_interests_placeholder` | manuscript-declarations | PENDING_AUTHOR_OR_EXTERNAL_ACTION | Provisional competing-interest sentence present in draft.tex, draft_siads_review.tex. | Confirm with authors and replace with final competing-interest declaration. |
| `funding_placeholder` | manuscript-declarations | PENDING_AUTHOR_OR_EXTERNAL_ACTION | Provisional funding sentence present in draft.tex, draft_siads_review.tex. | Confirm funding/no-funding statement with authors and insert final wording. |
| `author_metadata_incomplete` | front-matter | PENDING_AUTHOR_OR_EXTERNAL_ACTION | Current author block gives Yao Li affiliation/email only. | Confirm final author order, Jayleen affiliation/email, corresponding author, and ORCID choices. |
| `data_release_route_unfinalized` | data-availability | PENDING_AUTHOR_OR_EXTERNAL_ACTION | Data availability currently cites the GitHub repository but no immutable release tag or DOI. | Choose GitHub release-only or DOI-backed raw-data archive route; update data availability and rerun the gate. |
| `target_journal_confirmation` | journal-system | PENDING_AUTHOR_OR_EXTERNAL_ACTION | No target-journal choice can be proven from local manuscript text. | Authors must confirm target journal, article type, and whether SIADS review formatting is final. |
| `professional_similarity_screening` | external-service | PENDING_AUTHOR_OR_EXTERNAL_ACTION | No iThenticate/Turnitin/journal similarity report is available in the repository. | Run professional similarity/self-plagiarism screening after final edits. |
| `raw_data_doi_upload` | external-service | PENDING_AUTHOR_OR_EXTERNAL_ACTION | A local raw-data .tar.gz build exists, but no DOI/upload can be proven locally. | If the DOI route is chosen, upload the archive to Zenodo/OSF/journal storage and insert the DOI. |
| `final_post_edit_pdf_review` | journal-system | PENDING_AUTHOR_OR_EXTERNAL_ACTION | Current PDFs are local build artifacts before final author/journal edits. | After final declarations and release metadata are inserted, rerun the gate and visually inspect the final PDF. |

## Interpretation

This audit records author, journal, and external-service blockers.  A
`PENDING_AUTHOR_OR_EXTERNAL_ACTION` item is not a local numerical or
LaTeX failure, but the manuscript should not be formally submitted
until the item is resolved and the full local gate is rerun.
