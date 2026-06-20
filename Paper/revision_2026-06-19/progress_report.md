# Progress report — 2026-06-19

## Completed in this revision pass

1. Preserved the original extracted manuscript/code and created a revision
   working copy under `Paper/revision_2026-06-19/`.
2. Audited the draft and identified P0 model/data consistency issues.
3. Replaced the invalid June 18 flux workflow with a canonical
   Gibbs-preserving action-current simulator:
   `flux/NLS_flux_canonical.cpp`.
4. Added reproducible analysis scripts:
   - `flux/analyze_canonical_flux.py`
   - `flux/analyze_current_windows.py`
   - `flux/gibbs_mcmc_reference.py`
   - `flux/compare_gibbs_sde.py`
5. Ran validation gates:
   - warning-clean optimized C++ build;
   - sanitizer smoke test;
   - deterministic one-thread/two-thread reproducibility;
   - equal-temperature zero-current check;
   - hot/cold reversal check;
   - independent Gibbs MCMC profile check.
6. Ran production current-scaling experiments for `n=10,20,30,40` with
   corrected boundary noise and correct standard errors.
7. Updated `draft.tex` to:
   - correct the full-coordinate diffusion coefficient;
   - separate full phase variables from doubled bond-angle variables;
   - fix the Fourier-law sign convention;
   - use action-current terminology instead of heat-current terminology;
   - remove unsupported entropy/fluctuation-theorem claims;
   - replace the invalid `n^-1.75` claim with the corrected `n^-1.850`
     production result;
   - fix the stable phase-locking branch;
   - weaken unsupported “real spectrum/overdamped” claims;
   - insert the corrected scaling figure.
8. Restored local image resources needed by the short-chain figures and compiled
   the revision cleanly.
9. Added submission-critical manuscript sections:
   - English abstract;
   - conclusion and limitations;
   - data/code availability;
   - declarations for ethics, author contributions, competing interests,
     funding status, and AI-assisted preparation.
10. Replaced the remaining long-chain figure placeholders with reproducible
    generated figures:
    - `action_profiles.pdf/png`;
    - `cascade_embedding.pdf/png`;
    - `lte_residual_midchain.pdf/png`.
    The generation script is
    `Paper/revision_2026-06-19/scripts/generate_manuscript_figures.py`, and
    source metrics are recorded in
    `Paper/revision_2026-06-19/manuscript_figure_metrics.json`.
11. Replaced the phase-locking “derivation to be supplied” note with the
    reduced fixed-point quadratic
    `4(1+gamma^2) sin(theta)^2 + 4 gamma sin(theta) - 3 = 0` and the stable
    branch selection.
12. Recompiled the manuscript with TeX Live at that stage of the revision. The
    PDF then had 17 pages and the LaTeX log had no unresolved
    citations/references and no overfull or underfull box warnings. Later
    reproducibility-summary and submission-level additions increased the
    current build artifact to 20 pages.
13. Resolved two remaining rigor issues in the short-chain section:
    - removed the unsupported quantitative “about 15%” symmetry-breaking claim
      and made the figure caption explicitly qualitative;
    - added a source-traced eigen relaxation-rate sensitivity analysis in
      `scripts/analyze_eigen_fit_windows.py`, with results in
      `eigen_fit_sensitivity.json`, and replaced `eigenvalue_scatter.png`
      with a reproducible relaxation-rate diagnostic figure. The historical
      `lambda_R=-0.934` value is now described as an observable-dependent
      diagnostic rate, not a high-accuracy spectral-gap estimate.
14. Converted the reference list from an inline `thebibliography` block to
    `references.bib`, verified all eight references against publisher/arXiv/DOI
    metadata, updated `DobsonLiZhai` to its published CMS version, and recorded
    the audit in `integrity_audit_2026-06-19.md`.
15. Added `scripts/export_source_trace_metrics.py` and
    `source_trace_metrics.json`, which regenerate all LTE table values from
    histogram/profile files and export archived short-chain neural-network
    notebook outputs plus figure hashes. The draft now aligns the LTE
    `T_kin` values with this source-trace export and removes/downgrades
    short-chain quantitative claims whose notebook cells had no saved output
    (angular-width ratio, phase-locking peak table, middle-mode current
    balance).
16. Extracted short-chain neural-network diagnostics into
    `scripts/recompute_short_chain_nn_metrics.py` and reran them with
    `/usr/local/bin/python3` (TensorFlow 2.21.0 / Keras 3.14.0), generating
    `short_chain_nn_rerun_metrics.json`. The rerun reproduces the equilibrium
    validation errors and eigen-surrogate data-fit RMSE from the archived
    notebooks. It also confirms that the removed angular-width `≈2` claim is
    unsupported: the saved-model rerun gives `sigma3/sigma1≈1.03--1.13`, so the
    manuscript continues not to use that result as a claim.
17. Added `scripts/audit_manuscript_claims.py`, which builds
    `manuscript_claim_audit.json` and `manuscript_claim_audit.md` from the
    source-traced numerical artifacts. The current audit verifies 18/18 core
    numerical/data claims after the larger-chain and bath-temperature
    robustness updates and caught a notation issue in the LTE section:
    the equilibrium marginal is now written consistently with the manuscript's
    Gibbs convention as `exp[-H/(2T)]`.
18. Added a compact `Numerical reproducibility summary` table to the manuscript,
    mapping each numerical result family to its evidence bundle, automated or
    source-traced check, and retained claim scope. This addresses the previous
    recommendation to include a concise reproducibility table without adding
    new unsupported numerical claims.
19. Added `originality_spotcheck_2026-06-19.md`, a local web exact-phrase
    spot-check over short manuscript fragments and author/title queries. The
    spot-check found no visible exact external phrase reuse in returned
    titles/snippets and identified only expected topical sources, including the
    cited `HLNS` prior work. It remains a limited pre-screen, not a substitute
    for iThenticate/Turnitin-style professional checking.
20. Added `submission_readiness_checklist_2026-06-19.md` and appended current
    resolution-status updates to the initial audit, material inventory, and
    revision roadmap. The old documents now explicitly distinguish their
    initial-draft findings from the current revised-manuscript state.
21. Tightened the data/code availability statement so listed artifacts use
    repository-root paths where appropriate, then added
    `scripts/audit_availability_paths.py`, `availability_path_audit.json`, and
    `availability_path_audit.md`. The audit checks every `\path{...}` entry and
    every manuscript figure include; after the 2026-06-20 robustness and
    submission-readiness passes the current local result is 37/37 paths present
    and tracked. The LaTeX build
    remains clean with no unresolved
    citations/references and no overfull or underfull box warnings after the
    path-list formatting adjustment.
22. Added `scripts/build_submission_bundle_manifest.py` plus
    `submission_bundle_manifest.json` and `submission_bundle_manifest.md`. The
    manifest merges manuscript source, bibliography, figures, availability
    paths, claim-audit evidence, and handoff documents into one release
    checklist. It also records a raw-data limitation: 44 source-trace raw-data
    dependency records exist locally under untracked roots such as
    `Energy Cascade/`, `KDE/`, and `lte/`; these should be handled by a
    deliberate archive/Zenodo/OSF policy if the target journal requires full
    raw-data release.
23. Added `scripts/build_raw_data_archive_manifest.py`,
    `raw_data_archive_manifest.json`, and `raw_data_archive_manifest.md`. This
    converts the raw-data limitation into an actionable archive subset: 40
    unique source-trace raw files, all present locally, totaling 138,875,181
    bytes. The referenced subset is much smaller than the full local roots
    (`Energy Cascade/`, `KDE/`, and `lte/`), so a future DOI-backed supplement
    can be prepared deliberately without committing multi-GB directories to
    Git.
24. Added `scripts/run_submission_checks.py` plus
    `submission_checks_summary.json` and `submission_checks_summary.md`. Running
    the script with `--compile-latex` performs the local submission gate in one
    command: LaTeX compile/log scan, availability-path audit, manuscript-claim
    audit, raw-data archive manifest, and submission-bundle manifest. The
    current result is
    `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`.
25. Added `author_submission_action_packet_2026-06-19.md`, an author-facing
    packet that turns the remaining non-code blockers into fillable decisions:
    target journal/template, author metadata, funding, competing interests,
    CRediT-style contributions, professional similarity checking, raw-data DOI
    archive choice, and a generic cover-letter template.
26. Added `target_journal_shortlist_2026-06-19.md`, a preliminary shortlist
    based on official journal/publisher pages for SIADS, Physica D, Journal of
    Statistical Physics, Nonlinearity, Journal of Nonlinear Science, and Chaos.
    It records fit, risks, and template/data implications so the authors can
    choose a target before journal-specific conversion.
27. Added `submission_reproducibility_readme_2026-06-19.md`, a journal-neutral
    reviewer/editor entry point explaining the fast verification command, the
    purpose of each audit/manifest artifact, the raw-data archive convention,
    and the remaining author/journal-only items.
28. Added `pdf_layout_qa_2026-06-19.md`, a compiled-PDF layout QA record. The
    generic `article` PDF builds to 18 A4 pages with no scanned LaTeX warnings,
    overfull/underfull boxes, unresolved-reference markers, placeholder tokens,
    blank pages, obvious clipping, or visible figure/table rendering failures in
    the rendered-page inspection.
29. Added `scripts/audit_references.py` plus
    `reference_integrity_audit.json` and `reference_integrity_audit.md`. The
    audit verifies the local TeX/BibTeX citation graph: 8 BibTeX entries, 28
    citation commands, 29 citation uses, 8 unique cited keys, and zero dangling
    citations, orphan references, missing required fields, missing DOI/URL
    identifiers, missing recorded external verification sources, or duplicate
    BibTeX keys.
30. Added `scripts/build_submission_source_bundle.py` plus
    `submission_source_bundle_report.json` and
    `submission_source_bundle_report.md`. The script builds a source-only
    submission `.tar.gz` under `tmp/` from the release manifest, records an
    archive SHA-256 checksum, and deliberately excludes self-referential
    generated summaries and the large local raw-data roots. The current
    packaging run includes 259 regular files, has zero missing files, and
    records 44 local raw-data dependency records as
    excluded pending any DOI-backed raw-data archive decision.
31. Performed a submission-level strengthening pass on 2026-06-20:
    - added a Monte Carlo validation and uncertainty-protocol subsection to
      the manuscript;
    - added the LTE residual mesh figure generated from the
      `compare_residual.m` convention, including the requested `n=15`
      residual diagnostic;
    - added a timestep sensitivity table for the current estimator;
    - added finite-window current diagnostics for standardized current
      windows and $\tau\,\mathrm{Var}(\overline J_\tau)$;
    - reran the one-command submission gate, which now reports
      `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`
      with 37/37 path records available, 18/18 registered numerical claims
      verified, and 9 author/external submission items pending.
32. Added `siads_first_submission_packet_2026-06-20.md`, a target-specific
    preparation packet for the recommended first journal target.  It contains a
    SIADS-facing cover-letter draft, supplementary-material index, keywords and
    MSC candidates, data/code availability wording for two release routes, and
    a conversion checklist.  It is intentionally not a template conversion:
    SIADS still requires author confirmation before changing the manuscript
    class/style.
33. Added `draft_siads_review.tex`, a separate line-numbered SIADS
    review-preparation source.  The local TeX installation did not include the
    SIAM article class, so the file follows SIADS' non-SIAM-macro fallback with
    `lineno`, keywords, and MSC candidates.  It compiles cleanly to
    `tmp/paper_build/siads_review/draft_siads_review.pdf`. After the later
    robustness edits, the current compiled review PDF is 21 pages with SHA-256
    `1534e68e0617d0f44bc3d4b856d1ee608029ba3cc84c409a64802cfc93a4950b`.
34. Strengthened the manuscript narrative so the long-chain and short-chain
    results have a clearer division of labor: the long-chain simulations now
    explicitly provide macroscopic finite-size transport and LTE evidence,
    while the three-mode Fokker--Planck computations are framed as a
    low-dimensional microscope for stabilization, phase-locking, and
    qualitative slow-mode diagnostics.  The same wording was synchronized into
    `draft_siads_review.tex`, and both manuscript sources compile under the
    local LaTeX workflow.
35. Added a larger-chain current robustness check at `n=50` using the same
    canonical Gibbs-preserving current accumulator, `dt=5e-4`, measurement
    window `200`, and `1024` trajectories.  The run gives
    `E[J(50)] = 0.01851584685` with SE `0.00044158954`; adding this point to
    the four primary production lengths gives a diagnostic five-size exponent
    `-1.89449` with bootstrap 95% CI `[-1.91636,-1.87340]`.  A smaller
    fine-step pilot at `dt=2.5e-4` gives `E[J(50)] = 0.01879771710` with SE
    `0.00081439495`, only `1.52%` above the `dt=5e-4` run (`0.30` pooled
    standard errors), and an independent burn-in-10000 check gives
    `E[J(50)] = 0.01931242054` with SE `0.00085798987`.  The manuscript now
    records this as a robustness check rather than replacing the primary
    `n=10,20,30,40` exponent, because `n=50` is still a single larger-length
    extension and the fine-step run is a smaller pilot rather than a full
    production-resolution convergence study.

## Key validated numerical result

For the corrected Gibbs-preserving SDE at `T1=10`, `Tn=2`, `gamma=0.1`,
`dt=5e-4`, burn-ins `(1000,1280,2880,5120)`, and `1024` trajectories per chain
length:

| n | mean action current | SE |
|---:|---:|---:|
| 10 | 0.3925219606 | 0.0018902080 |
| 20 | 0.1191693526 | 0.0009195305 |
| 30 | 0.0545731139 | 0.0006191849 |
| 40 | 0.0297475540 | 0.0004827205 |

Power-law fit:

- `E[J(n)] = 28.7457 n^-1.85008`
- `R^2 = 0.998013`
- bootstrap 95% exponent CI: `[-1.87034, -1.83049]`

## Remaining blockers before journal submission

- Funding statement still needs author confirmation. The current manuscript
  explicitly says funding information was not supplied in the available
  materials.
- Author-contribution and competing-interest declarations should be confirmed
  by the authors before submission.
- Target journal, article type, and author metadata should be filled in using
  `author_submission_action_packet_2026-06-19.md`; candidate venues and
  conversion implications are summarized in
  `target_journal_shortlist_2026-06-19.md`.
- Professional plagiarism/self-plagiarism screening is still required before
  formal submission. The local short-fragment web spot-check is clean within its
  sampled scope, but it is not a corpus-scale similarity report.
- Rerun the final citation/data/claim/path audit after author declarations and
  any last journal-format edits are settled. The current local core numerical
  audit passes 18/18 checks, and the current availability-path audit passes
  37/37 checks. The preferred one-command local gate is now
  `python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex`.
  The reviewer/editor-facing navigation file is
  `submission_reproducibility_readme_2026-06-19.md`; the compiled-PDF layout
  QA snapshot is `pdf_layout_qa_2026-06-19.md`; the local citation/reference
  graph audit is `reference_integrity_audit.md`; the source-only submission
  packaging report is `submission_source_bundle_report.md`.
- Decide whether the final release should include only the audited derived
  artifacts in GitHub or also a DOI-backed raw-data archive. The bundle
  manifest currently flags 44 local raw-data dependency records that are not
  git-tracked; the raw-data archive manifest deduplicates them to 40 unique
  files totaling 138,875,181 bytes.
- Optional but recommended before final release: package a full TensorFlow
  retraining environment if the target journal expects re-training
  reproducibility beyond saved-model inference.
- Optional numerical strengthening: if time permits before submission, turn the
  `n=50` robustness check into a full larger-length convergence study by
  increasing the fine-step run to production resolution and, if
  computationally feasible, adding one further length such as `n=60`.
