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
| Short-chain saved-model diagnostics | Rerun from saved TensorFlow/Keras models | `scripts/recompute_short_chain_nn_metrics.py`; `short_chain_nn_rerun_metrics.json` |
| Data/code and figure path availability | Passing local path audit | `scripts/audit_availability_paths.py`; `availability_path_audit.json`; `availability_path_audit.md` |
| Submission/release bundle | Passing for tracked release files; raw-data limitation recorded | `scripts/build_submission_bundle_manifest.py`; `submission_bundle_manifest.json`; `submission_bundle_manifest.md` |
| Source-only submission archive | Passing packaging dry run; archive written under `tmp/` | `scripts/build_submission_source_bundle.py`; `submission_source_bundle_report.json`; `submission_source_bundle_report.md` |
| Minimal raw-data archive plan | Passing local raw-file manifest; archive not yet uploaded | `scripts/build_raw_data_archive_manifest.py`; `raw_data_archive_manifest.json`; `raw_data_archive_manifest.md` |
| One-command local gate | Passing with raw-data archive limitation | `scripts/run_submission_checks.py`; `submission_checks_summary.json`; `submission_checks_summary.md` |
| Reproducibility entry point | Prepared for reviewer/editor navigation | `submission_reproducibility_readme_2026-06-19.md` |
| Compiled-PDF layout QA | Generic article PDF checked locally | `pdf_layout_qa_2026-06-19.md` |
| Author/journal action packet | Prepared; requires author completion | `author_submission_action_packet_2026-06-19.md` |
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
as an extrapolation of the transport exponent.  The latest one-command gate
reports `PASS_WITH_LOCAL_RAW_DATA_LIMITATION` with 34/34 availability path
records present and git-tracked, 14/14 registered numerical claims verified,
and zero missing release-bundle files.

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
2. Confirm `submission_checks_summary.md` reports
   `PASS_WITH_LOCAL_RAW_DATA_LIMITATION` or better.
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
   or raw-data dependency changes; if a Zenodo/OSF upload is made, verify the
   archived file set against `raw_data_archive_manifest.md`.
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
  - `submission_bundle_manifest.md`
  - `submission_source_bundle_report.md`
  - `raw_data_archive_manifest.md`
  - `submission_checks_summary.md`
  - `submission_reproducibility_readme_2026-06-19.md`
  - `pdf_layout_qa_2026-06-19.md`
  - `author_submission_action_packet_2026-06-19.md`
  - `target_journal_shortlist_2026-06-19.md`
  - `siads_first_submission_packet_2026-06-20.md`
  - `submission_readiness_checklist_2026-06-19.md`
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
