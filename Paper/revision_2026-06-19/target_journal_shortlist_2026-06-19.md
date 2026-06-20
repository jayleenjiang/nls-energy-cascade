# Target-journal shortlist — updated 2026-06-20

Scope: preliminary target-journal shortlist for
`Paper/revision_2026-06-19/draft.tex`, based on official journal or publisher
pages checked on 2026-06-19 and spot-verified again for the top candidates on
2026-06-20.  This is not a submission decision; it is a decision aid for the
authors.

Current manuscript profile:

- numerical/theoretical study of a stochastic energy-cascade chain derived from
  a resonant NLS toy model;
- methods/results mix: stochastic dynamics, nonequilibrium steady states,
  long-chain Monte Carlo, local-equilibrium diagnostics, short-chain
  Fokker--Planck/neural diagnostics, reproducibility audits;
- current local gate status:
  `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`;
- remaining non-local decisions: author declarations, professional
  similarity/self-plagiarism check, journal template, and possible DOI-backed
  raw-data archive.

## Shortlist summary

| Candidate | Fit | Main reason to consider | Main risk / preparation needed |
|---|---:|---|---|
| SIAM Journal on Applied Dynamical Systems (SIADS) | High | Strong fit for mathematical analysis/modeling of dynamical systems and computational reproducibility. | Needs SIAM-style cover letter/PDF submission, inline figures, and clear supplementary-material index if raw/auxiliary files are submitted. |
| Physica D: Nonlinear Phenomena | High | Strong fit for nonlinear phenomena, Hamiltonian/integrable systems, wave motion, turbulence, and data-driven dynamical systems. | Elsevier-specific formatting/data statement should be checked in the submission system before conversion. |
| Journal of Statistical Physics | Medium-high | Good fit for nonequilibrium/statistical-mechanics framing and stochastic steady-state/current scaling. | Need Springer Nature formatting/declarations; paper should keep statistical-mechanics contribution prominent. |
| Nonlinearity | Medium-high | Good fit for nonlinear dynamics and mathematical/numerical nonlinear science. | IOP supplementary-file size limits make a repository/DOI archive important if raw data are included. |
| Journal of Nonlinear Science | Medium | Broad nonlinear-science venue that explicitly allows numerical simulations when they give broader insight. | Higher bar for broad multidisciplinary framing; introduction may need more non-specialist motivation. |
| Chaos: An Interdisciplinary Journal of Nonlinear Science | Medium | Strong nonlinear-science audience and interest in comprehensible manifestations of nonlinear phenomena. | Manuscript would need more interdisciplinary/plain-language framing; AIP formatting requirements should be checked. |

## Candidate notes

### 1. SIAM Journal on Applied Dynamical Systems (SIADS)

Official pages checked:

- Journal page: <https://www.siam.org/publications/siam-journals/siam-journal-on-applied-dynamical-systems/>
- Instructions for authors: <https://epubs.siam.org/journal/siads/instructions-for-authors>

Why it fits:

- The official SIADS page describes the journal as publishing mathematical
  analysis and modeling of dynamical systems and applications across physical,
  engineering, life, and social sciences.
- The manuscript is centered on a dynamical system, stochastic forcing,
  numerical steady states, and reproducibility of computational diagnostics.

Preparation implications:

- SIADS asks authors to submit a manuscript and cover letter in PDF format via
  the submission system, and figures should be embedded inline in the
  manuscript.
- SIADS has explicit supplementary-material handling; if raw data or code are
  submitted as supplement, prepare an index describing each item and why it is
  included.
- The current generic `article` draft is acceptable for internal review, but a
  SIAM conversion should use the SIAM article style and check abstract/reference
  conventions before submission.

Verdict: best first target if the authors want an applied-math/dynamical-systems
venue and are willing to keep the reproducibility/supplementary package tidy.

### 2. Physica D: Nonlinear Phenomena

Official pages checked:

- Elsevier journal page: <https://www.elsevier.com/locate/physd>
- Elsevier journal listing/scope page:
  <https://shop.elsevier.com/journals/physica-d-nonlinear-phenomena/0167-2789>

Why it fits:

- Elsevier describes Physica D as publishing theoretical and experimental work,
  techniques, and ideas that advance understanding of nonlinear phenomena.
- The listed scope includes nonlinear-system methods, wave motion,
  hydrodynamics/turbulence, integrable and Hamiltonian systems, and
  data-driven dynamical systems, all close to the manuscript's NLS-cascade and
  numerical-NESS frame.

Preparation implications:

- Before conversion, verify the current Elsevier submission-system requirements
  for data availability, declaration statements, and preferred LaTeX format.
- The manuscript's current framing should emphasize nonlinear phenomena,
  Hamiltonian/wave-turbulence motivation, and finite-size numerical discovery.

Verdict: strong physics/nonlinear-phenomena target, especially if the authors
prefer a nonlinear-science audience over a purely applied-math audience.

### 3. Journal of Statistical Physics

Official pages checked:

- Submission guidelines:
  <https://link.springer.com/journal/10955/submission-guidelines>
- Springer Nature LaTeX support:
  <https://www.springernature.com/gp/authors/campaigns/latex-author-support>

Why it fits:

- The paper studies nonequilibrium steady states, stochastic baths,
  finite-size current scaling, local-equilibrium diagnostics, and limitations
  around fluctuation/entropy-production claims.
- Springer's guidelines for this journal support LaTeX submissions and require
  editable source files plus compiled PDF.

Preparation implications:

- Convert to Springer Nature LaTeX template if selected.
- Confirm competing-interest, funding, data/code availability, and declaration
  statements using Springer wording.
- Strengthen the statistical-mechanics narrative in the introduction/abstract
  if the target is chosen.

Verdict: good target if the authors want a statistical-mechanics framing, but
the paper should avoid sounding like a purely numerical methods report.

### 4. Nonlinearity

Official page checked:

- IOP Publishing Support page:
  <https://publishingsupport.iopscience.iop.org/journals/nonlinearity/>

Why it fits:

- The manuscript is a nonlinear stochastic/dynamical systems study derived from
  dispersive PDE resonance structure.
- IOP states that common TeX/LaTeX variants are acceptable and that use of the
  IOP class file is helpful but not essential.

Preparation implications:

- IOP notes that supplementary material/data files can be limited by file size
  and recommends a data repository for larger data sets.  The current minimal
  raw-data archive subset is about 139 MB, so a DOI-backed data repository may
  be safer than direct journal supplement upload.
- If selected, check whether numerical evidence and finite-size limitations are
  framed strongly enough for a nonlinear-analysis readership.

Verdict: plausible target if the authors emphasize nonlinear dynamics and
mathematical structure, with repository-backed data handling.

### 5. Journal of Nonlinear Science

Official page checked:

- Submission guidelines:
  <https://link.springer.com/journal/332/submission-guidelines>

Why it fits:

- The journal explicitly includes theory, experimentation, algorithms,
  numerical simulations, and applications, provided the work is creative and
  sound.
- Its author guidance asks computational papers to validate results and provide
  sufficient information for reproduction; this aligns with the current audit
  package.

Preparation implications:

- The journal emphasizes broad nonlinear-science significance.  The paper's
  introduction may need an extra broad-context pass if selected.
- The guidelines say submissions should be in LaTeX or PDF and ask authors to
  suggest suitable editorial-board members.

Verdict: possible but more ambitious; use if the authors can sharpen the broad
conceptual contribution beyond the NLS-cascade niche.

### 6. Chaos: An Interdisciplinary Journal of Nonlinear Science

Official pages checked:

- AIP journal page: <https://pubs.aip.org/aip/cha>
- AIP author instructions: <https://publishing.aip.org/resources/researchers/author-instructions/>
- AIP Chaos editorial policies page:
  <https://pubs.aip.org/aip/cha/pages/policies>

Why it fits:

- Chaos is explicitly devoted to nonlinear phenomena and aims for explanations
  comprehensible to researchers from a broad range of disciplines.
- Recent issue examples include turbulent and stochastic/nonlinear dynamics
  topics, indicating broad numerical-dynamics coverage.

Preparation implications:

- The manuscript would need an accessibility pass: less specialist notation in
  the introduction, more explanatory framing for the cascade and NESS physics,
  and perhaps a stronger graphical abstract or conceptual figure if requested.
- AIP accepts Word or LaTeX in its general author instructions; verify
  journal-specific Chaos requirements before conversion.

Verdict: plausible interdisciplinary venue if the authors want a broader
nonlinear-science readership and are willing to adapt the exposition.

## Practical recommendation

For the current manuscript with minimal additional rewriting:

1. **First practical target:** SIADS.
2. **Strong physics/nonlinear alternative:** Physica D.
3. **Statistical-mechanics alternative:** Journal of Statistical Physics.

If the authors choose **Journal of Nonlinear Science** or **Chaos**, plan a
separate broad-audience introduction/framing pass before template conversion.

## 2026-06-20 official-policy spot check

The recommendation above remains unchanged after checking the current official
author pages for the leading candidates:

- **SIADS** remains the cleanest first target.  SIAM asks for manuscript and
  cover letter PDFs at submission, requires inline figures, encourages SIAM
  macros, and generally expects manuscripts not to exceed 40 pages or 10 MB
  without cover-letter justification.  SIADS also requires an index if
  supplementary materials are submitted.  The current 20-page generic article
  draft is therefore comfortably within the length/file-size envelope; the main
  remaining SIADS-specific work is template/line-number conversion and a
  supplementary-material index if the authors choose to submit raw data or code
  as SIADS supplement.
- **Physica D** remains the strongest nonlinear-physics alternative.  Elsevier
  asks for editable source files, encourages the Elsevier LaTeX template, and
  requires a generative-AI disclosure when AI tools were used in manuscript
  preparation.  The current manuscript already has an AI-assisted-preparation
  statement, but it should be adapted to Elsevier's section title and wording
  if Physica D is selected.
- **Journal of Statistical Physics** is a good third target if the paper is
  framed around NESS/current scaling/local equilibrium.  Springer Nature asks
  for LaTeX source plus compiled PDF and requires a data availability statement
  for original research articles.  The current data/code availability section
  is close, but a DOI-backed raw-data archive would make this target safer.
- **Nonlinearity** remains plausible but the data-package logistics are less
  convenient for this project: IOP accepts common TeX/LaTeX variants, but its
  supplementary files are limited to 50 MB each and 150 MB combined, with larger
  data recommended for a repository.  The current minimal raw-data subset is
  about 139 MB, so a repository/DOI route is preferable if Nonlinearity is
  chosen.

**Actionable recommendation:** unless the authors have a strong physics-journal
preference, prepare the next version for SIADS first.  Do not convert the class
file until the authors confirm this target; the current generic article PDF is
better for internal review and coauthor comments.

## Template/conversion checklist after target selection

1. Confirm journal and article type.
2. Confirm whether the journal requires anonymized review.
3. Convert LaTeX class/style only after target selection.
4. Insert final author affiliations, ORCID IDs, corresponding author, funding,
   competing interests, and contribution statements.
5. Decide data-release route:
   - GitHub-only derived reproducibility artifacts; or
   - GitHub plus DOI-backed minimal raw-data archive from
     `raw_data_archive_manifest.md`.
6. Rerun:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

7. Run professional originality/self-plagiarism checking on the final PDF/source.
