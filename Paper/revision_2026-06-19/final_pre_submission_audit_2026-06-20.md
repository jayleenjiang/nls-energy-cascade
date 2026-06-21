# Final local pre-submission audit — 2026-06-20

Scope: `Paper/revision_2026-06-19/draft.tex` and the synchronized
`draft_siads_review.tex` after the larger-chain `n=50`, production-size
`n=60`, bath-temperature and thermostat-coupling robustness extensions,
fit-window sensitivity, and fine-timestep checks.

This audit separates items that are locally complete from items that still
require author, journal, or external-service action.  It is intentionally
conservative: the paper should not be formally submitted until the
author-required items below are resolved in the manuscript source and the local
gate is rerun.

## Local evidence status

Latest local command:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Latest result: `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`.

Key gate numbers:

| Gate | Status |
|---|---|
| LaTeX/log scans | PASS for `draft.tex` and `draft_siads_review.tex`, `issues=0` in both logs |
| Compiled PDF artifact audit | PASS, `2/2` local PDFs verified with page counts, sizes, and SHA-256 checksums |
| Availability/path audit | PASS, `38/38` paths present, `0` untracked required files |
| Numerical claim audit | PASS, `22/22` claims verified |
| Reference integrity audit | PASS, `8` cited BibTeX entries, `0` dangling citations |
| Author/submission-field audit | `AUTHOR_CONFIRMATION_PENDING`, `9` pending author/external items |
| Raw-data archive manifest | PASS locally, `40` unique referenced files, `138,875,181` bytes |
| Raw-data archive build | PASS locally, upload-ready `.tar.gz` prepared under `tmp/`, DOI/upload still external |
| Submission bundle manifest | `PASS_WITH_LOCAL_RAW_DATA_LIMITATION`, `121` release files, `0` missing, `0` untracked release files |
| Submission metadata consistency audit | PASS, `21` handoff metadata checks over `9` documents |
| Source-only bundle dry run | PASS, `317` included files |
| SIADS cover-letter template build | PASS locally; template PDF contains placeholders and is not final |
| Journal upload package build | PASS locally, SIADS repository-route package written under `tmp/journal_upload_package/runs/` |

The SIADS review-preparation source is now included directly in the
one-command local gate:

- source: `Paper/revision_2026-06-19/draft_siads_review.tex`
- PDF: `tmp/paper_build/siads_review/draft_siads_review.pdf`
- local compile status: PASS

## Reviewer-style local stress test

### Core thesis

The core numerical thesis is now appropriately scoped.  The main transport
claim remains the production fit over `n=10,20,30,40`; the new `n=50` and
`n=60` experiments are explicitly described as larger-chain robustness checks,
and the `T1=8,Tn=4` production-resolution run plus the `gamma=0.05` and
`gamma=0.2` production-resolution runs are described as finite-size parameter
robustness checks.  Together with the fit-window sensitivity analysis and a
production-resolution `n=50` fine-step check at `dt=2.5e-4`, this avoids the most likely
reviewer objections that larger lengths or alternate bath parameters have been
over-promoted to an asymptotic law.

### Current-scaling vulnerability

The remaining current-scaling vulnerability is not a local inconsistency: it is
the usual finite-size limitation.  The manuscript now says this directly.  A
reviewer could still ask for production-size fine-step data at `n=60`; that
would strengthen the paper but is not needed to defend the present finite-size
claim as written, and the `n=50` fine-step request has now been addressed.
The added fit-window sensitivity table shows
that the existing `n=10`--`60` windows do not drift toward the Fourier exponent
over the simulated range.  The production-size `n=60` extension gives a
diagnostic six-length exponent near `-1.93`, so it is promoted to manuscript
robustness evidence while the primary exponent remains the original
four-length production fit.

A production-resolution bath-parameter robustness check has also been recorded
under `experiments/flux_validation/parameter_robustness_2026-06-20/`.  The
`T1=8,Tn=4` run gives exponent `-1.75098` with bootstrap 95% CI
`[-1.77964,-1.72269]` and maximum stationarity statistic `1.73684` paired
standard errors.  It supports the faster-than-Fourier finite-size trend under a
second bath-temperature pair, while the manuscript explicitly stops short of a
systematic parameter sweep.

A production-resolution thermostat-coupling robustness check has also been
recorded under
`experiments/flux_validation/gamma_robustness_2026-06-21/`.  At
`T1=10,Tn=2`, the `gamma=0.05` run gives exponent `-1.65035` with bootstrap
95% CI `[-1.66794,-1.63333]`, while the `gamma=0.2` run gives exponent
`-1.99149` with bootstrap 95% CI `[-2.01710,-1.96682]`.  The maximum
split-window stationarity statistics are `1.14405` and `1.74247` paired
standard errors.  These checks support the manuscript's narrower statement
that the faster-than-Fourier finite-size trend is not tied only to
`gamma=0.1`, while still stopping short of a systematic two-parameter sweep.

### Short-chain mechanism vulnerability

The short-chain Fokker--Planck/eigenfunction material is now framed as a
mechanistic microscope rather than proof of the long-chain exponent.  The
neural solver claims are limited to saved-model inference and source-traced
diagnostics.  The added short-chain solver-diagnostics table consolidates the
Gibbs-slice, angular-asymmetry, phase-locking, current-balance, and
backward-Monte-Carlo checks while explicitly reserving high-accuracy spectral
and quantitative NESS claims.  This is the safer journal posture.

### Local-equilibrium vulnerability

The LTE section distinguishes pair-marginal agreement from strict local Gibbs
structure, and the residual mesh diagnostic is included.  The claim audit
source-traces the table values and convention.  No local blocker remains.

### Entropy-production and large-deviation vulnerability

The manuscript explicitly does not claim bath-energy entropy production,
Gallavotti--Cohen symmetry, or asymptotic current large deviations.  This
removes a likely overclaim.

## Author-required manuscript edits before formal submission

The following are the remaining text-level blockers inside the current
manuscript.  They require author confirmation and therefore have not been
silently converted into final declarations.

### 1. Author approval/contribution sentence

Current text in `draft.tex`:

```tex
Both authors should review and approve the final submitted version.
```

Replace after confirmation with a final factual statement, for example:

```tex
Both authors reviewed and approved the final submitted version.
```

or with the exact target-journal authorship wording.

### 2. Competing-interest declaration

Current text:

```tex
No competing interests are declared in the materials supplied for this draft.
```

Replace after author confirmation with either:

```tex
The authors declare no competing interests.
```

or a complete disclosure.

### 3. Funding declaration

Current text:

```tex
Funding information was not supplied in the current manuscript materials and
should be completed by the authors before submission if required by the target
journal.
```

Replace after author confirmation with either:

```tex
The authors received no external funding for this work.
```

or the complete grant/funding statement required by the target journal.

### 4. Author metadata

Confirm final author order, affiliations, corresponding author, email, and
ORCID identifiers if requested by the target journal.  The current source
contains only Yao Li's affiliation/email in a footnote.

## External or journal-system items still required

1. Select final target journal and article type.
2. Run professional plagiarism/self-plagiarism screening, such as
   iThenticate/Turnitin or the target journal's required equivalent.
3. Decide whether the GitHub release is sufficient or whether the 40-file
   raw-data subset should be uploaded to Zenodo/OSF with a DOI.
4. If submitting to SIADS, decide whether to keep the line-numbered fallback
   source for review or convert to the SIAM macro package after installing the
   journal class locally.
5. Rerun the full local gate after inserting final declarations and any
   journal-specific formatting changes.

## Local conclusion

No remaining blocker is a local numerical, citation, path, source-bundle, or
LaTeX compilation failure.  The remaining blockers are author-confirmed
declaration wording, journal-format decisions, professional similarity
screening, and the optional DOI-backed raw-data archive.
