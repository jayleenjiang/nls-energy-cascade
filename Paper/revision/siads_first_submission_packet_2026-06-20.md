# SIADS-first submission packet — 2026-06-20

Scope: practical preparation packet if the authors choose **SIAM Journal on
Applied Dynamical Systems (SIADS)** as the first target for
`Paper/revision_2026-06-19/draft.tex`.

This packet does **not** mean the manuscript has been converted or submitted.
It translates the current generic, locally checked manuscript into the concrete
SIADS-facing tasks that remain after author confirmation.

Official SIADS pages checked on 2026-06-20 and refreshed against current
SIAM/SIADS policy pages on 2026-06-22.  The detailed policy refresh is recorded
in `target_journal_policy_refresh_2026-06-22.md`.

- Journal page:
  <https://www.siam.org/publications/siam-journals/siam-journal-on-applied-dynamical-systems/>
- Instructions for authors:
  <https://epubs.siam.org/journal/siads/instructions-for-authors>
- Editorial policy:
  <https://epubs.siam.org/journal/siads/editorial-policy>
- SIAM AI policy:
  <https://epubs.siam.org/artificial-intelligence>

## Current readiness snapshot

- Recommended first target: SIADS.
- Current manuscript format: generic `article`; the current
  SIADS review-preparation PDF is 28 pages.
- SIADS review-preparation source:
  `Paper/revision_2026-06-19/draft_siads_review.tex`.
- SIADS review-preparation PDF build:
  `tmp/paper_build/siads_review/draft_siads_review.pdf`.
- SIADS review-preparation PDF SHA-256:
  `f85aa910afdadfb2112589191405a3d75f5cac978cd3ac6add87284149be4bfc`.
- Current local gate:
  `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`.
- Current manuscript figures are embedded inline.
- Current code-verifiable numerical claim audit: 23/23 claims verified.
- Current availability-path audit: 43/43 paths present and required files
  git-tracked.
- Remaining non-local blockers: target confirmation, author metadata,
  funding/competing-interest declarations, professional similarity check, and
  raw-data archive decision.
- Author/submission-field audit:
  `AUTHOR_CONFIRMATION_PENDING` with 9 pending author/external items.

## SIADS-specific checklist

Complete these only after the authors confirm SIADS as the target.

1. **Manuscript format**
   - Use `draft_siads_review.tex` as the current SIADS review-preparation
     source.  It is a line-numbered copy of `draft.tex` with keywords and MSC
     codes added.
   - If SIAM macros are available before submission, optionally convert the
     source to SIAM/SIADS style with review mode.  The local TeX installation
     used for this pass did not include SIAM's article class, so the prepared
     version follows SIADS' non-SIAM-macro line-numbering fallback.
   - Preserve inline figures.
   - Recompile and visually inspect the converted PDF.

2. **Length and file size**
   - Keep the manuscript under the SIADS general 40-page / 10 MB expectation.
   - If the converted PDF exceeds either threshold, justify it in the cover
     letter or move nonessential material to supplement.

3. **Front matter**
   - Add final affiliations, corresponding author, and emails.
   - Add ORCID IDs if the authors want them displayed.
   - Add keywords and MSC codes.
   - Recheck that the abstract remains one paragraph and under 250 words after
     any target-specific edits.

4. **Declarations**
   - Confirm funding statement.
   - Confirm competing-interest statement.
   - Confirm author-contribution statement.
   - Confirm AI-assisted-preparation wording against SIAM's May 2026 AI policy,
     including the authors' responsibility for all content.

5. **Supplementary material**
   - Decide whether to submit supplementary files to SIADS or instead cite a
     GitHub/Zenodo/OSF archive.
   - If supplementary files are submitted with SIADS, include an index listing
     each item, its description, and why it should accompany the paper.

6. **Final checks before upload**
   - Run:

     ```sh
     python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
     ```

   - Run professional similarity/self-plagiarism screening.
   - Confirm the manuscript is not under consideration elsewhere.
   - Upload manuscript PDF and cover letter PDF through the SIADS submission
     system.

## Suggested SIADS keywords and MSC codes

Author should confirm these before submission.

Suggested keywords:

- resonant nonlinear Schrödinger equation
- energy cascade
- nonequilibrium steady state
- stochastic Hamiltonian system
- anomalous transport
- local equilibrium
- Fokker--Planck equation

Suggested MSC candidates:

- 35Q55 — NLS-like nonlinear dispersive equations
- 37M05 — computational methods in dynamical systems
- 60H10 — stochastic ordinary differential equations
- 82C31 — stochastic methods in nonequilibrium statistical mechanics
- 82C70 — transport processes

## SIADS cover letter draft

Replace bracketed fields before use.
The same text is also available as a compilable LaTeX template in
`siads_cover_letter_template.tex`; its local template-PDF build is audited in
`siads_cover_letter_template_build.md`.  The template PDF is not final and
must not be submitted until all bracketed fields are replaced and authors
approve the cover letter.

```text
Dear Editors,

We are pleased to submit our manuscript entitled
"Numerical study of an energy cascade model derived from a dispersive
equation" for consideration in SIAM Journal on Applied Dynamical Systems.

The manuscript studies a stochastic energy-cascade chain derived from the
resonant toy model for the cubic nonlinear Schrödinger equation.  The rigorous
three-mode nonequilibrium steady-state theory for this model leaves longer
chains open; our work provides a reproducible numerical study of those
long-chain nonequilibrium steady states and returns to the three-mode system
for detailed Fokker--Planck diagnostics.

The main contributions are:

1. formulation and validation of a Gibbs-preserving two-bath stochastic model;
2. high-throughput Monte Carlo evidence that the finite-chain action current
   satisfies E[J(n)] ≈ 28.75 n^{-1.850} over n = 10, 20, 30, 40 under the
   reported parameters, with n = 50 and production-size n = 60 robustness
   checks plus bath-temperature and thermostat-coupling robustness runs;
3. local-equilibrium diagnostics showing that long-chain pair marginals are
   close to rescaled equilibrium marginals at the local kinetic temperature,
   while strict local Gibbs structure fails;
4. source-traced three-mode neural Fokker--Planck diagnostics, including
   equilibrium validation, qualitative symmetry breaking, phase locking, and a
   carefully scoped slow-mode diagnostic;
5. a reproducibility package with scripted claim audits, path audits,
   release-bundle manifests, and raw-data archive manifests.

The paper is well suited to SIADS because it combines stochastic/deterministic
dynamical-systems modeling, numerical invariant-measure computation, and
reproducible diagnostics for a nonlinear energy-cascade model derived from a
dispersive equation.

The manuscript is not under consideration elsewhere.  All authors have
approved the submitted version.  [Insert funding statement, competing-interest
statement, and any SIADS-required AI-use disclosure.]

[If supplementary material is submitted: We also submit supplementary material
consisting of [brief list].  An index of the supplementary materials is
included.]

Sincerely,

[Corresponding Author Name]
[Affiliation]
[Email]
```

## Proposed SIADS supplementary-material index

Use this only if the authors decide to submit supplementary material directly
with SIADS.  If a DOI-backed archive is used instead, replace this table with
the archive DOI and a shorter statement.

| Item | Proposed filename / location | Description | Justification |
|---|---|---|---|
| S1 | `submission_reproducibility_readme_2026-06-19.md` | Reviewer/editor navigation file for the reproducibility package. | Gives referees a concise map of the audit scripts, generated artifacts, and local gate command. |
| S2 | `submission_source_bundle_report.md` and source-only archive generated by `scripts/build_submission_source_bundle.py` | Manifest and source-only archive of manuscript source, figures, scripts, validation artifacts, and checksums. | Provides a compact source package for reproducing the manuscript build and verifying included files. |
| S3 | `experiments/flux_validation/` | Current-scaling validation artifacts, production summaries, larger-chain robustness, bath-parameter robustness, current-window diagnostics, and validation report. | Supports the main action-current scaling, robustness checks, and finite-window current diagnostics. |
| S4 | `manuscript_claim_audit.md`, `availability_path_audit.md`, `reference_integrity_audit.md` | Machine-readable and human-readable audit outputs. | Documents that registered numerical claims, local paths, and references are internally consistent. |
| S5 | `raw_data_archive_manifest.md` plus DOI-backed raw-data archive, if created | Minimal source-traced raw-data subset: 42 files, 151,605,557 bytes, paths preserved under `raw_data/`. | Needed only if SIADS/referees require raw histograms/model files beyond the tracked derived artifacts. |

## Data/code availability wording for SIADS draft

Use after the authors choose the release route.

**GitHub-only route:**

> The manuscript source, analysis scripts, generated figures, validation
> artifacts, and audit manifests are available in the project GitHub repository
> at [repository release URL].  The repository includes scripts for checking
> manuscript paths, references, registered numerical claims, and the submission
> source bundle.  Large local raw-data roots are summarized in
> `raw_data_archive_manifest.md` and can be archived separately if requested.

**GitHub + DOI-backed raw-data route:**

> The manuscript source, analysis scripts, generated figures, validation
> artifacts, and audit manifests are available in the project GitHub repository
> at [repository release URL].  The raw files needed to reproduce the
> source-traced numerical tables and figures are archived at [Zenodo/OSF DOI],
> with paths and SHA-256 checksums listed in `raw_data_archive_manifest.md`.

## Do not submit until these are filled

- Target journal confirmed as SIADS by both authors.
- Final author order, affiliations, emails, and corresponding author confirmed.
- Funding and competing-interest statements confirmed.
- Author-contribution statement confirmed.
- Professional similarity/self-plagiarism check completed.
- Raw-data release route chosen.
- SIAM/SIADS format conversion completed and final gate rerun.
