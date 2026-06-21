# Journal submission readiness checklist — 2026-06-19

Scope: final local handoff checklist for
`Paper/revision_2026-06-19/draft.tex`.  This file records the remaining
non-code items that must be supplied or confirmed by the authors before the
paper can be called journal-ready.

## Current local status

| Area | Status | Evidence |
|---|---|---|
| Main manuscript | Locally revised and source-traced | `draft.tex`; `progress_report.md`; `integrity_audit_2026-06-19.md` |
| Framework/source copy | Preserved | `paper_draft_1.tex`; `material_inventory.md`; original archive/backups |
| References | Verified and converted to BibTeX | `references.bib`; `integrity_audit_2026-06-19.md` |
| Citation/reference integrity | Passing local structural audit | `scripts/audit_references.py`; `reference_integrity_audit.json`; `reference_integrity_audit.md` |
| Core numerical claims | Passing local claim audit | `scripts/audit_manuscript_claims.py`; `manuscript_claim_audit.json`; `manuscript_claim_audit.md` |
| Figure provenance | Source-traced for manuscript-generated figures; LTE mesh figure now included | `scripts/generate_manuscript_figures.py`; `manuscript_figure_metrics.json`; `report_assets/compare_residual_mesh.pdf` |
| Monte Carlo validation and finite-window diagnostics | Added to manuscript and checked by local gate | `experiments/flux_validation/validation_report.md`; `experiments/flux_validation/production_dt5e-4/current_windows_window_statistics.csv`; `draft.tex` |
| Gamma-robustness production path | Smoke-tested only; not manuscript evidence | `scripts/run_gamma_robustness_smoke.py`; `gamma_robustness_smoke_report.md` |
| Short-chain saved-model diagnostics | Rerun from saved TensorFlow/Keras models | `scripts/recompute_short_chain_nn_metrics.py`; `short_chain_nn_rerun_metrics.json` |
| Data/code and figure path availability | Passing local path audit | `scripts/audit_availability_paths.py`; `availability_path_audit.json`; `availability_path_audit.md` |
| Compiled-PDF artifact metadata | Passing for generic and SIADS local PDFs | `scripts/audit_compiled_pdfs.py`; `compiled_pdf_artifact_audit.json`; `compiled_pdf_artifact_audit.md` |
| Submission/release bundle | Passing for tracked release files; raw-data limitation recorded | `scripts/build_submission_bundle_manifest.py`; `submission_bundle_manifest.json`; `submission_bundle_manifest.md` |
| Source-only submission archive | Passing packaging dry run; archive written under `tmp/` | `scripts/build_submission_source_bundle.py`; `submission_source_bundle_report.json`; `submission_source_bundle_report.md` |
| Journal upload package builder | Prepared; local package written under `tmp/` by the one-command gate | `scripts/build_journal_upload_package.py`; `journal_upload_file_index_2026-06-20.md` |
| SIADS cover-letter template | Prepared and locally compiled; not final until bracketed fields are replaced and authors approve | `siads_cover_letter_template.tex`; `scripts/build_siads_cover_letter_template.py`; `siads_cover_letter_template_build.md` |
| Submission metadata consistency | Passing local audit for handoff PDF hashes, page counts, and bundle counts | `scripts/audit_submission_metadata_consistency.py`; `submission_metadata_consistency_audit.json`; `submission_metadata_consistency_audit.md` |
| Minimal raw-data archive plan | Passing local raw-file manifest; upload-ready local archive build prepared but not uploaded | `scripts/build_raw_data_archive_manifest.py`; `scripts/build_raw_data_archive.py`; `raw_data_archive_manifest.json`; `raw_data_archive_manifest.md`; `raw_data_archive_build_report.md` |
| One-command local gate | Passing locally with explicit author-confirmation and raw-data archive limitations; compiles both the generic manuscript and SIADS review source | `scripts/run_submission_checks.py`; `submission_checks_summary.json`; `submission_checks_summary.md`; `author_submission_fields_audit.md` |
| Reproducibility entry point | Prepared for reviewer/editor navigation | `submission_reproducibility_readme_2026-06-19.md` |
| Compiled-PDF layout QA | Generic article PDF checked locally | `pdf_layout_qa_2026-06-19.md` |
| Author/journal action packet | Prepared; requires author completion | `author_submission_action_packet_2026-06-19.md` |
| Final author/journal information request | Prepared as a concise fillable form mapped to `author_submission_fields_template.json` | `final_author_submission_fields_request_2026-06-20.md` |
| Author-field application workflow | Prepared; dry-run-first and backup-on-apply | `author_submission_fields_template.json`; `scripts/apply_author_submission_fields.py` |
| Target-journal shortlist | Prepared from official pages; requires author choice | `target_journal_shortlist_2026-06-19.md` |
| SIADS-first submission packet | Prepared for the recommended first target; requires author confirmation before use | `siads_first_submission_packet_2026-06-20.md` |
| Originality pre-screen | Clean within sampled web-query scope | `originality_spotcheck_2026-06-19.md` |

## 2026-06-20 submission-level update

The current manuscript now includes a dedicated Monte Carlo validation and
uncertainty-protocol subsection, the `n=15,25,50` LTE residual mesh diagnostic
from the `compare_residual.m` convention, a timestep sensitivity table for the
current estimator, and finite-window current diagnostics.  A subsequent
manuscript pass clarified the division of labor between the long-chain
transport/LTE evidence and the three-mode Fokker--Planck mechanism study, so
the short-chain section is framed as a low-dimensional microscope rather than
as an extrapolation of the transport exponent.  A further current-scaling
robustness pass added a fit-window sensitivity table around the `n=50` and
`n=60` larger-chain runs.  A subsequent bath-parameter robustness pass added a
production-resolution `T1=8,Tn=4` current-scaling check over
`n=10,20,30,40`.  The latest one-command gate reports
`PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION` with
37/37 availability path records present and git-tracked, 19/19 registered
numerical claims verified, 9 author/external submission items pending, and zero
missing release-bundle files.

## 2026-06-20 larger-chain robustness update

A new `n=50` current run was added under
`experiments/flux_validation/larger_n_pilot_2026-06-20/` using the canonical
current accumulator, `dt=5e-4`, burn-in `8000`, measurement window `200`, and
`1024` trajectories.  It gives `E[J(50)] = 0.01851584685` with SE
`0.00044158954`, and the diagnostic fit over `n=10,20,30,40,50` gives exponent
`-1.89449` with bootstrap 95% CI `[-1.91636,-1.87340]`.  A smaller
fine-timestep pilot at `dt=2.5e-4` gives `E[J(50)] = 0.01879771710` with SE
`0.00081439495`, differing from the `dt=5e-4` result by `1.52%` or `0.30`
pooled standard errors.  The fit-window sensitivity analysis gives exponents
`-1.85008` on `n=10,20,30,40`, `-1.89449` on `n=10,20,30,40,50`,
`-1.92956` on `n=10,20,30,40,50,60`, and `-2.05926` on the tail
`n=20,30,40,50,60`.  The manuscript uses these as robustness checks only, not
as the primary quoted exponent, because a systematic larger-length and
fine-timestep convergence study is still outside the revision scope.

Subsequent `n=60` pilots were recorded under
`experiments/flux_validation/larger_n60_pilot_2026-06-20/`.  The initial
128-trajectory pilot gives `E[J(60)] = 0.01189574358` with SE `0.00122139546`;
the medium 256-trajectory pilot gives `E[J(60)] = 0.01361149053` with SE
`0.00081058642`.  A production-size 1024-trajectory run gives
`E[J(60)] = 0.01244829643` with SE `0.00041661977` and stationarity statistic
`-0.524` paired SE.  Including the production-size `n=60` run in a diagnostic
six-length fit gives exponent `-1.92956` with bootstrap 95% CI
`[-1.95424,-1.90603]`.  This is now suitable as manuscript robustness
evidence, while the primary exponent remains the four-length production fit.

## 2026-06-20 bath-parameter robustness update

A production-resolution bath-temperature robustness check was added under
`experiments/flux_validation/parameter_robustness_2026-06-20/`.  The clean
pilot set `T1=8,Tn=4` was upgraded to `1024` trajectories per length over
`n=10,20,30,40`, with the same `dt=5e-4`, measurement window `200`, and
burn-in schedule as the primary production runs.  The fitted exponent is
`-1.75098` with bootstrap 95% CI `[-1.77964,-1.72269]` and log-fit
`R^2=0.99844`; the maximum split-window stationarity statistic is `1.73684`
paired standard errors.  The manuscript uses this as a robustness check that
the faster-than-Fourier finite-size decay is not tied only to the single
temperature pair `T1=10,Tn=2`, not as a systematic parameter sweep.

## Author confirmations still required

These items cannot be completed from the local code/data alone.

Use `author_submission_action_packet_2026-06-19.md` as the fillable handoff
document for these confirmations and `target_journal_shortlist_2026-06-19.md`
as the preliminary venue-decision aid.

1. **Target journal and template**
   - Confirm target journal.
   - Decide whether to keep generic `article` format for first circulation or
     convert to a journal-specific LaTeX class.
   - Check journal-specific requirements for declarations, data availability,
     AI disclosure, reference style, figure resolution, and supplementary
     material.

2. **Author metadata**
   - Confirm final author order.
   - Confirm affiliations for all authors.
   - Confirm corresponding author and email.
   - Add ORCID identifiers if the target journal requests them.

3. **Author-contribution statement**
   - Confirm whether the current statement is accurate:
     Jayleen Jiang performed numerical experiments, assembled computational
     artifacts, and drafted the manuscript; Yao Li supervised the project and
     contributed to model formulation, theoretical framing, and interpretation.
   - If the journal uses CRediT taxonomy, map contributions to those categories.

4. **Funding and competing interests**
   - Supply funding award names/numbers or state explicitly that there was no
     external funding.
   - Confirm competing-interest statement.

5. **External originality check**
   - Run a professional plagiarism/self-plagiarism check such as
     iThenticate/Turnitin after final author and journal-format edits.
   - Compare against prior drafts, reports, preprints, theses, and any
     submitted manuscripts by the author team.

6. **Data/code release decision**
   - Decide whether the GitHub branch is sufficient for review or whether to
     create an immutable release/tag.
   - If the journal requires archival data, prepare a Zenodo/OSF release and
     replace the GitHub-only availability statement with DOI-backed language.
   - Decide whether large local raw-data roots (`Energy Cascade/`, `KDE/`, and
     `lte/`) should be archived outside GitHub. The current bundle manifest
     records 44 local source-trace raw-data dependency records that exist
     locally but are not git-tracked; the raw-data archive manifest
     deduplicates these to 40 unique files totaling 138,875,181 bytes.

7. **Neural-network reproducibility level**
   - Decide whether saved-model inference reproducibility is enough.
   - If the journal requests full retraining, archive the TensorFlow/Keras
     training environment and notebooks/scripts needed to regenerate the saved
     models from raw training data.

## Final local gates to rerun after author edits

Run these only after the author-supplied items above have been inserted.

1. Run
   `python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex`
   from the repository root.
2. Confirm `submission_checks_summary.md` reports no local numerical, path,
   reference, bundle, or LaTeX failures. Before author-only items are resolved
   the expected status is
   `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`; after
   those items are resolved, rerun until no author-confirmation placeholders
   remain.
3. Compile `draft.tex` from a clean build directory if the runner was not run
   with `--compile-latex`.
4. Check the LaTeX log for unresolved references/citations and overfull or
   underfull box warnings.
5. Visually inspect every page of the compiled PDF.
6. Rerun `scripts/audit_manuscript_claims.py` and confirm all registered
   numerical/data claims pass.
7. Rerun `scripts/audit_availability_paths.py` and confirm every manuscript
   `\path{...}` entry and figure include exists in the release branch or
   archive.
8. Rerun `scripts/build_submission_bundle_manifest.py` and confirm the tracked
   release bundle has no missing or untracked required files. If a raw-data
   archive is created, rerun the manifest/checklist against that archive.
9. Rerun `scripts/build_submission_source_bundle.py` and confirm
   `submission_source_bundle_report.md` reports `PASS`.
10. Rerun `scripts/build_raw_data_archive_manifest.py` if any source-trace JSON
    or raw-data dependency changes; rerun
    `scripts/build_raw_data_archive.py` to produce the local upload `.tar.gz`;
    if a Zenodo/OSF upload is made, verify the archived file set against
    `raw_data_archive_manifest.md` and `raw_data_archive_build_report.md`.
11. Recheck `references.bib` for dangling or orphan citation keys.
12. Update `integrity_audit_2026-06-19.md` or create a new dated final audit
   snapshot with the final compile/audit results.

## Suggested final submission bundle

- Manuscript source: `draft.tex`
- Bibliography: `references.bib`
- Compiled manuscript PDF
- Figure files referenced by `draft.tex`
- Reproducibility/readme material:
  - `progress_report.md`
  - `integrity_audit_2026-06-19.md`
  - `reference_integrity_audit.md`
  - `manuscript_claim_audit.md`
  - `availability_path_audit.md`
  - `compiled_pdf_artifact_audit.md`
  - `submission_bundle_manifest.md`
  - `submission_source_bundle_report.md`
  - `siads_cover_letter_template_build.md`
  - `submission_metadata_consistency_audit.md`
  - `raw_data_archive_manifest.md`
  - `raw_data_archive_build_report.md`
  - `author_submission_fields_audit.md`
  - `submission_checks_summary.md`
  - `submission_reproducibility_readme_2026-06-19.md`
  - `pdf_layout_qa_2026-06-19.md`
  - `author_submission_action_packet_2026-06-19.md`
  - `final_author_submission_fields_request_2026-06-20.md`
  - `author_submission_fields_template.json`
  - `target_journal_shortlist_2026-06-19.md`
  - `siads_first_submission_packet_2026-06-20.md`
  - `journal_upload_file_index_2026-06-20.md`
  - `final_submission_decision_sheet_2026-06-20.md`
  - `submission_readiness_checklist_2026-06-19.md`
  - `final_pre_submission_audit_2026-06-20.md`
- Optional supplementary archive:
  - current-scaling validation artifacts under `experiments/flux_validation/`
  - figure and source-trace metrics JSON files
  - short-chain saved-model rerun script/output
  - eigen-fit sensitivity script/output

## Current blocker summary

The paper is locally much closer to submission quality than the initial draft:
major numerical/model inconsistencies have been corrected, unsupported claims
have been removed or demoted, and core code-verifiable numerical claims pass
the local audit.  The remaining blockers are not new simulations; they are
author/journal/external-verification items: author metadata and declarations,
professional originality screening, target-journal formatting, and a final
post-edit audit.
