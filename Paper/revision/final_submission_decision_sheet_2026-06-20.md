# Final submission decision sheet — 2026-06-20

Scope: practical author/advisor decision sheet for the current
`Paper/revision_2026-06-19/draft.tex` manuscript and the synchronized
`draft_siads_review.tex` review-preparation source.

This sheet separates what is locally complete from what still requires an
author, journal, or external-service decision.  It should be used before
changing the manuscript from a locally verified draft into a formally submitted
paper.

## One-line status

The manuscript is locally reproducible and internally audited, but it is not
formally submission-ready until author declarations, target-journal choices,
similarity screening, and the data-release route are confirmed.

## Current local evidence

Latest local gate:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Latest status: `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`.

Key numbers from the current gate:

| Item | Current value |
|---|---:|
| LaTeX/log checks | `draft.tex` PASS; `draft_siads_review.tex` PASS |
| Availability paths checked | 43/43 |
| Registered numerical claims verified | 23/23 |
| Cited BibTeX entries | 11 |
| Dangling citations | 0 |
| Compiled PDF artifact audit | 2/2 PDFs verified |
| Release-bundle files | 130 |
| Missing release files | 0 |
| Source-only bundle included files | 326 |
| Minimal raw-data subset | 42 files, 151,605,557 bytes |
| Author/submission-field audit | 9 pending author/external items |

The SIADS review-preparation source compiles locally:

- source: `Paper/revision_2026-06-19/draft_siads_review.tex`
- PDF: `tmp/paper_build/siads_review/draft_siads_review.pdf`
- PDF SHA-256:
  `f85aa910afdadfb2112589191405a3d75f5cac978cd3ac6add87284149be4bfc`
- PDF pages: 28

The current upload-facing local file index is
`journal_upload_file_index_2026-06-20.md`.  It lists the generic manuscript
PDF, SIADS review PDF, source-only archive, and minimal raw-data archive with
local paths and SHA-256 checksums.  The current local PDF metadata is also
recorded in `compiled_pdf_artifact_audit.md`.

## Scientific readiness decision

Local scientific/reproducibility status: **yes, locally ready for advisor
review as a near-submission manuscript**.

The current manuscript now has:

1. a corrected Gibbs-preserving long-chain stochastic model;
2. production current scaling over `n=10,20,30,40`;
3. larger-chain robustness at `n=50` and production-size `n=60`;
4. timestep and burn-in/stationarity diagnostics;
5. a second bath-temperature robustness run at `T1=8,Tn=4`;
6. thermostat-coupling robustness runs at `gamma=0.05` and `gamma=0.2`;
7. LTE pair-marginal diagnostics with the residual mesh figure;
8. carefully scoped finite-window current statistics;
9. short-chain Fokker--Planck/neural/eigenfunction diagnostics framed as a
   mechanism study rather than proof of the long-chain exponent;
10. reproducibility, path, claim, reference, bundle, and raw-data manifests.

Formal submission status: **not yet**.  The remaining blockers are not local
numerical failures; they are author/journal/external-verification items.

## Author decisions required before submission

### 1. Target journal and article type

Recommended first practical target from the current packet: SIADS.

Files to use:

- `target_journal_shortlist_2026-06-19.md`
- `target_journal_policy_refresh_2026-06-22.md`
- `siads_first_submission_packet_2026-06-20.md`

Decision needed:

- submit first to SIADS, or choose another target;
- if SIADS is chosen, decide whether to keep the current line-numbered review
  source or convert to official SIAM multimedia macros before upload.

### 2. Author metadata

Confirm:

- final author order;
- affiliations for all authors;
- corresponding author and email;
- ORCID identifiers if desired or required.

### 3. Author-contribution wording

Current draft wording still contains a pre-submission placeholder:

```tex
Both authors should review and approve the final submitted version.
```

Replace only after confirmation, for example:

```tex
Both authors reviewed and approved the final submitted version.
```

or with the exact target-journal authorship wording.

### 4. Funding and competing interests

Current draft wording is deliberately conservative because the information was
not supplied locally.

Confirm and replace with final factual statements, for example:

```tex
The authors declare no competing interests.
```

and either:

```tex
The authors received no external funding for this work.
```

or the complete grant/funding statement required by the target journal.

Implementation aid: use
`final_author_submission_fields_request_2026-06-20.md` as the concise
author-facing form.  After authors confirm these fields, copy
`author_submission_fields_template.json` to `author_submission_fields.json`,
fill the fields, run `scripts/apply_author_submission_fields.py` in dry-run
mode, and use `--apply` only after it reports `ready_to_apply=true`.

### 5. Similarity/self-plagiarism screening

Run a professional screen such as iThenticate, Turnitin, or the target
journal's required equivalent after final author and journal-format edits.

The local originality spot-check is useful but does not replace professional
screening.

### 6. Data-release route

Choose one route before replacing the data-availability placeholder:

| Route | When to use | Manuscript consequence |
|---|---|---|
| GitHub release only | fastest review route if the journal accepts repository artifacts | cite the immutable repository release URL |
| GitHub + Zenodo/OSF DOI | strongest archival route and recommended for journal submission | cite both repository release and raw-data DOI |
| Local raw-data available on request | weakest route; use only if archive is impossible | keep the raw-data limitation explicit |

The strongest route is GitHub release plus a DOI-backed raw-data archive using
the 42-file, 151,605,557-byte subset in `raw_data_archive_manifest.md`.  A
local upload-ready `.tar.gz` can now be built with
`scripts/build_raw_data_archive.py`; the latest local build is recorded in
`raw_data_archive_build_report.md`.  This still does not create a DOI or upload
the archive.

## Optional additional studies

These are not blockers for the current finite-size manuscript as written, but
they would strengthen a more ambitious version:

1. production or medium-production larger-chain current runs at `n=70,80,100`;
2. a broader systematic two-parameter bath sweep beyond the checked
   `T1=10,Tn=2`, `T1=8,Tn=4`, `gamma=0.05`, `gamma=0.1`, and `gamma=0.2`
   finite-size runs;
3. a fuller LTE residual-norm convergence table beyond the current
   source-traced mesh-slice RMS audit;
4. a separate bath-energy accumulator for entropy production and
   Gallavotti--Cohen diagnostics.

Do not add these results to the manuscript unless they are run with the same
source-tracing and audit standard used for the current validation artifacts.

## Final local gate after author edits

After inserting final declarations, author metadata, target-specific formatting,
and any release/DOI information, rerun:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Submission should proceed only if the gate remains
`PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION` or
better before author-only items are resolved, and then a fully author-confirmed
passing status after those items are inserted. The compiled PDF should also be
visually inspected page by page.
