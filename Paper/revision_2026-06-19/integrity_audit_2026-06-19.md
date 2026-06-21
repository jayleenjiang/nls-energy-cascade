# Integrity audit snapshot — 2026-06-19

Scope: `Paper/revision_2026-06-19/draft.tex` after the BibTeX refactor and
short-chain diagnostic tightening.

Verdict: **PASS WITH BLOCKING NOTES for submission**.  Bibliographic existence,
ghost citations, main current-scaling claims, LTE table values, and the
manuscript-used short-chain neural-network diagnostics pass the checks below.
The local short-fragment originality spot-check found no visible exact external
phrase reuse in returned web snippets, but it is not a substitute for a
professional similarity report.  Before journal submission, the authors must
still confirm funding / author-contribution declarations and run the final
external originality check.

## Reference verification

All eight references in `references.bib` were checked against publisher,
arXiv, PMLR, DOI, or institutional metadata pages.  The hand-written
`thebibliography` block was replaced with BibTeX so future formatting changes
can be handled by journal style files.

| Key | Verdict | Verification source | Notes |
|---|---|---|---|
| `CKSTT` | VERIFIED | Springer article page, DOI `10.1007/s00222-010-0242-2`: <https://link.springer.com/article/10.1007/s00222-010-0242-2> | Confirms authors, title, Invent. Math. 181, 39--113 (2010). |
| `HLNS` | VERIFIED | arXiv record: <https://arxiv.org/abs/2505.16018> | Corrected title to arXiv's singular form, `Non-equilibrium steady state for a three-mode energy cascade model`; added arXiv DOI. |
| `ZhaiDobsonLi` | VERIFIED | PMLR page: <https://proceedings.mlr.press/v145/zhai22a.html> | Confirms PMLR 145:568--597 (2022) and official BibTeX metadata. |
| `Li2019` | VERIFIED | International Press / DOI metadata and PDF record, DOI `10.4310/CMS.2019.v17.n4.a9`: <https://intlpress.com/JDetail/1806262739393794050> | Confirms Commun. Math. Sci. 17(4):1045--1059 (2019). |
| `DobsonLiZhai` | VERIFIED | International Press page, DOI `10.4310/CMS.2022.v20.n3.a8`: <https://link.intlpress.com/JDetail/1806261569648545793> | Updated from arXiv-only citation to published Commun. Math. Sci. 20(3):803--827 (2022). |
| `GallavottiCohen` | VERIFIED | APS DOI page: <https://link.aps.org/doi/10.1103/PhysRevLett.74.2694> | Confirms Phys. Rev. Lett. 74(14):2694--2697 (1995). |
| `LepriLiviPoliti` | VERIFIED | DOI/institutional metadata: <https://abdn.elsevierpure.com/en/publications/thermal-conduction-in-classical-low-dimensional-lattices/> | Confirms Physics Reports 377(1):1--80 (2003), DOI `10.1016/S0370-1573(02)00558-6`. |
| `Nazarenko` | VERIFIED | Springer book page: <https://link.springer.com/book/10.1007/978-3-642-15942-8> | Confirms Lecture Notes in Physics 825, Springer, 2011, DOI and ISBN. |

## Ghost citation check

Command:

```sh
python3 - <<'PY'
import re
from pathlib import Path
tex=Path('Paper/revision_2026-06-19/draft.tex').read_text()
bib=Path('Paper/revision_2026-06-19/references.bib').read_text()
tex_nc='\n'.join(line.split('%',1)[0] for line in tex.splitlines())
keys=[]
for m in re.finditer(r'\\cite\w*\s*(?:\[[^\]]*\]\s*){0,2}\{([^}]*)\}', tex_nc):
    keys.extend(k.strip() for k in m.group(1).split(',') if k.strip())
bibkeys=re.findall(r'@\w+\s*\{\s*([^,]+),', bib)
print('dangling citations', sorted(set(keys)-set(bibkeys)))
print('orphan bib entries', sorted(set(bibkeys)-set(keys)))
PY
```

Result:

- dangling citations: `[]`
- orphan BibTeX entries: `[]`

## Claim/data traceability snapshot

| Claim family | Manuscript locations | Evidence inspected | Status |
|---|---:|---|---|
| Corrected action-current scaling `E[J(n)] = 28.7457 n^-1.85008`, `R^2 = 0.998013`, bootstrap CI `[-1.87034,-1.83049]` | abstract, current-scaling section, conclusion | `experiments/flux_validation/production_manifest.md`; `experiments/flux_validation/validation_report.md`; `production_dt5e-4/flux_primary_scaling.json` | VERIFIED for finite-size numerical claim over `n=10,20,30,40`. |
| Gibbs-preserving flux implementation and validation gates | current-scaling section | frozen source/binary hashes in `production_manifest.md`; equal-temperature and hot/cold reversal gates in `validation_report.md` | VERIFIED for reported production workflow. |
| Long-chain action profiles and sample counts | action-profile section and figure caption | `manuscript_figure_metrics.json`; script `scripts/generate_manuscript_figures.py` | VERIFIED for `n=25,50,100` plotted profile values and sample counts. |
| Cascade embedding figure | CKSTT embedding paragraph and figure | `manuscript_figure_metrics.json`; source CSVs under `Energy Cascade/Rect Plots/` | VERIFIED as a generated visualization from local CSVs. |
| Mid-chain LTE residual figure | `fig:lte-resid` | `manuscript_figure_metrics.json`; script `scripts/generate_manuscript_figures.py`; hist files `lte/n50 data/simd_n50_j24.hist` and `lte/n100 data/n100_j48.hist` | VERIFIED for displayed residual slices and recorded fit metrics. |
| LTE residual mesh diagnostic | `fig:lte-resid-mesh` | `report_assets/compare_residual_mesh.pdf`; local `compare_residual_mesh_metrics.txt`; source convention `compare_residual.m` | VERIFIED as a structural residual diagnostic for `n=15,25,50`; used qualitatively, not as a replacement for the fixed weighted-core estimator in `tab:lte`. |
| Full LTE table (`tab:lte`) slopes/R²/temperature comparisons | LTE section | `source_trace_metrics.json`; script `scripts/export_source_trace_metrics.py`; hist/profile source hashes recorded for each row. | VERIFIED for all table rows after aligning `T_kin` entries to the exported pair-averaged profile values. |
| Equilibrium Fokker--Planck validation errors `1.5%, 2.6%, 6.4%`; short-chain density figures; manuscript-used short-chain NN diagnostics | short-chain Fokker--Planck section | `source_trace_metrics.json`; `short_chain_nn_rerun_metrics.json`; `scripts/export_source_trace_metrics.py`; `scripts/recompute_short_chain_nn_metrics.py`; archived notebook outputs from `KDE/4:15_NN/FKE_5d_NLS.ipynb`; figure SHA256 records for `eq_validation.png`, `neq_density.png`, `symmetry_breaking.png`, and `Q1_slices.png`. | VERIFIED for archived equilibrium-error values, figure provenance, and saved-model TensorFlow rerun (`/usr/local/bin/python3`, TensorFlow 2.21.0, Keras 3.14.0). The rerun reproduces the equilibrium errors and eigen surrogate RMSE, supports qualitative symmetry breaking, and confirms the removed angular-width `≈2` claim is unsupported (`sigma3/sigma1≈1.03--1.13`). Middle-mode current balance has an 8.9% residual in the rerun and is retained only as a solver diagnostic, not a manuscript claim. |
| Eigen relaxation diagnostic `lambda_diag≈-0.934` and window sensitivity | eigenfunction section | `eigen_fit_sensitivity.json`; script `scripts/analyze_eigen_fit_windows.py`; source `KDE/NLS_backward_Y_train.txt` | VERIFIED as observable-dependent diagnostic rate, not a spectral-gap claim. |
| Symmetry breaking | stabilization section | `short_chain_nn_rerun_metrics.json`; manuscript now uses qualitative wording; old unsupported `15%` claim removed. | VERIFIED as a qualitative diagnostic. The saved-model rerun reports masked and unmasked asymmetry metrics, but the manuscript intentionally avoids a standalone percentage because the value is sensitive to low-density regions and network normalization. |
| Core numerical manuscript claims | abstract, claim-evidence map, long-chain section, LTE table, LTE residual decomposition, short-chain section, short-chain solver-diagnostics appendix, reproducibility summary, data availability | `manuscript_claim_audit.json`; `manuscript_claim_audit.md`; script `scripts/audit_manuscript_claims.py` | VERIFIED for 22/22 code-verifiable numerical/data claims after the larger-chain, bath-temperature, thermostat-coupling, LTE residual-decomposition, short-chain solver-diagnostics, and introductory claim-evidence-map updates. This audit also caught and corrected the LTE notation so the local equilibrium marginal uses the manuscript convention `exp[-H/(2T)]`, not the ambiguous `exp(-H/beta)`. |
| Citation/reference integrity | local TeX citation graph, BibTeX completeness, DOI/URL presence, and recorded external source URLs | `reference_integrity_audit.json`; `reference_integrity_audit.md`; script `scripts/audit_references.py` | VERIFIED structurally for 11/11 BibTeX entries and 11/11 unique cited keys: zero dangling citations, orphan references, missing required fields, missing DOI/URL identifiers, missing recorded external sources, or duplicate BibTeX keys. External source URLs were checked manually on 2026-06-19 and 2026-06-21 and recorded in the audit. |
| Data/code and figure path availability | data/code availability statement and all manuscript figure includes | `availability_path_audit.json`; `availability_path_audit.md`; script `scripts/audit_availability_paths.py` | VERIFIED locally for 38/38 referenced paths after the larger-chain, bath-parameter, thermostat-coupling, LTE-mesh, and finite-window current updates. This check proves local path existence and hashes, but should be rerun against the final release branch/archive. |
| Compiled PDF artifacts | generic manuscript PDF and SIADS review-preparation PDF under `tmp/` | `compiled_pdf_artifact_audit.json`; `compiled_pdf_artifact_audit.md`; script `scripts/audit_compiled_pdfs.py` | VERIFIED locally for 2/2 compiled PDFs: page counts, byte sizes, SHA-256 checksums, and TeX log byte summaries match the current local build artifacts. |
| Submission/release bundle manifest | manuscript source, bibliography, figures, availability paths, claim-audit evidence, and handoff documents | `submission_bundle_manifest.json`; `submission_bundle_manifest.md`; script `scripts/build_submission_bundle_manifest.py` | PASS WITH LOCAL RAW-DATA LIMITATION: the tracked release bundle has no missing required files after staging, but 44 local raw-data dependency records referenced by source-trace JSON exist only in local roots such as `Energy Cascade/`, `KDE/`, and `lte/`. A DOI-backed raw-data archive is still needed if the journal requires full raw-data release. |
| Source-only submission archive | package the tracked release files and tracked validation directory files into a checksumed `.tar.gz` under `tmp/` | `submission_source_bundle_report.json`; `submission_source_bundle_report.md`; script `scripts/build_submission_source_bundle.py` | PASS: the packaging dry run creates a source-only archive with zero missing files. The archive intentionally excludes self-referential generated summaries and the large local raw-data roots; raw data remain governed by `raw_data_archive_manifest.md` and a future DOI-backed archive decision. |
| Minimal raw-data archive manifest | source-trace raw files under `Energy Cascade/`, `KDE/`, and `lte/` | `raw_data_archive_manifest.json`; `raw_data_archive_manifest.md`; script `scripts/build_raw_data_archive_manifest.py` | VERIFIED locally for 40/40 unique referenced raw files, total 138,875,181 bytes. This is a preparation manifest only: it does not upload an archive, but it identifies a compact raw-data subset much smaller than the full local raw roots. |
| One-command local submission checks | LaTeX compile/log scans for `draft.tex` and `draft_siads_review.tex`, compiled-PDF artifact audit, availability audit, claim audit, reference audit, raw-data manifest, release bundle manifest, and source-bundle packaging | `submission_checks_summary.json`; `submission_checks_summary.md`; script `scripts/run_submission_checks.py` | PASS WITH AUTHOR CONFIRMATION PENDING AND LOCAL RAW-DATA LIMITATION when run with `--compile-latex`: both LaTeX logs have zero flagged issues, the PDF artifact audit passes 2/2 PDFs, availability audit passes 38/38, claim audit passes 22/22, reference audit passes, raw-data archive manifest passes 40/40, source-bundle packaging passes, and the tracked release bundle has no missing or untracked required files. |
| Reproducibility entry point | reviewer/editor navigation for the submission package | `submission_reproducibility_readme_2026-06-19.md` | PREPARED: the fast verification path, artifact map, raw-data archive convention, and non-local limitations are summarized in one journal-neutral entry point. |
| Compiled-PDF layout QA | generic `article` PDF layout and updated figure-page rendering checks | `pdf_layout_qa_2026-06-19.md` | PASS WITH LIMITATIONS: the 26-page A4 PDF has no scanned LaTeX-log warnings, and the newly added LTE mesh/current-diagnostic pages render without visible clipping after white-background conversion. A target-journal proof pass is still required after template conversion. |
| Author/journal action packet | target journal, author metadata, declarations, raw-data DOI decision, professional originality check, and cover letter | `author_submission_action_packet_2026-06-19.md` | PREPARED: remaining non-code blockers have been converted into a fillable author-facing packet. Completion still requires author input and external checks. |
| Target-journal shortlist | preliminary venue selection based on official journal/publisher pages | `target_journal_shortlist_2026-06-19.md` | PREPARED: candidate venues and conversion implications are documented. Completion still requires author selection of a target journal and any journal-specific template conversion. |
| Originality spot-check | sampled abstract, methods/results, short-chain, conclusion, and author/title queries | `originality_spotcheck_2026-06-19.md`; exact-phrase web queries over short manuscript fragments | PASS WITH LIMITATION: no exact external phrase reuse was visible in returned titles/snippets. Expected public `HLNS` prior work was found and is already cited. Professional plagiarism/self-plagiarism screening is still required before formal submission. |

## AI research failure mode checklist

| Mode | Status | Evidence / note |
|---|---|---|
| 1. Implementation bug passing self-review | CLEAR for canonical flux experiment, LTE table export, and saved-model short-chain NN inference; QUALIFIED for full NN retraining | Flux code has frozen hashes, sanitizer smoke run, deterministic threading check, equal-temperature and hot/cold gates. LTE table values are regenerated from hist/profile files. Short-chain NN values are exported from archived notebook outputs and rerun from saved Keras models in `scripts/recompute_short_chain_nn_metrics.py`; model retraining itself is not rerun here. |
| 2. Hallucinated citation | CLEAR | 8/8 references verified; no dangling/orphan citation keys. |
| 3. Hallucinated experimental result | CLEAR for flux/eigen/action-profile/LTE-table/LTE-residual and manuscript-used short-chain NN claims | Main scaling, gamma robustness, eigen diagnostics, action profiles, LTE table entries, LTE residual decomposition, and the introductory claim-evidence map are tied to JSON/manifest outputs. Short-chain NN equilibrium-error, qualitative symmetry, phase-locking branch, eigen-surrogate fit values, and the new solver-diagnostics table are source-traced through archived outputs and the saved-model rerun. The consolidated claim-audit script verifies 22/22 core numerical/data claims against source artifacts. The middle-mode current balance is source-traced only as a solver diagnostic, not as a manuscript claim. Unsupported unarchived quantitative claims were removed or downgraded. |
| 4. Shortcut reliance | CLEAR / not applicable | Numerical SDE simulation; no ML benchmark shortcut claim. Neural FP solver is now described diagnostically and validated against Gibbs where possible. |
| 5. Bug reframed as insight | CLEAR for corrected current-scaling narrative | Unsupported GC/entropy and spectral-gap claims were demoted; the paper labels finite-size and diagnostic limitations explicitly. |
| 6. Methodology fabrication | CLEAR for flux, LTE, and saved-model short-chain NN workflows; QUALIFIED for full NN retraining | Flux methods match manifest. LTE table computation is scripted. Short-chain NN workflow is traced to notebooks, model/data hashes, figure hashes, and a clean non-notebook TensorFlow rerun script. Full retraining of the neural networks is not part of the current release audit. |
| 7. Frame-lock | CLEAR WITH LIMITATION | Current framing is explicitly numerical/finite-size and no longer claims asymptotic theorem or full GC test. Remaining limitations are acknowledged. |

## Remaining required actions before calling the paper journal-ready

1. Confirm funding information and author-contribution / competing-interest
   declarations with the authors.
2. Run a professional plagiarism/self-plagiarism check outside this repository
   before formal submission; the local short-snippet web spot-check is clean
   within its sampled scope but is not a replacement for a professional
   similarity report.
3. Rerun the final citation/data/claim/path audit after author declarations and
   any last journal-format edits are finalized. The current local core
   numerical claim audit passes 22/22 checks and the availability-path audit
   passes 38/38 checks; the one-command runner
   `scripts/run_submission_checks.py --compile-latex` should be repeated after
   final edits.
4. Decide the raw-data release policy. The local source-trace artifacts point
   to large untracked raw-data roots (`Energy Cascade/`, `KDE/`, and `lte/`);
   however, the minimal raw-data archive manifest identifies 40 unique
   referenced files totaling 138,875,181 bytes. If the journal requires raw
   data rather than audited derived artifacts, create an archival DOI-backed
   supplement from that manifest and rerun the bundle manifest against it.
5. Optional but recommended: archive a full neural-network retraining
   environment if the journal requires re-training rather than saved-model
   inference reproducibility.
