# Author and journal submission action packet — 2026-06-19

Scope: author-facing packet for finishing
`Paper/revision_2026-06-19/draft.tex`.  The local numerical and path checks are
now scripted, but the items below require author or target-journal decisions.

Current local gate command:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Current local gate status: `PASS_WITH_LOCAL_RAW_DATA_LIMITATION`.

This means the manuscript compiles locally, the LaTeX log has no flagged
overfull/underfull/reference/citation warnings, registered numerical claims
pass the local audit, and release-bundle files are present and git-tracked.  It
does **not** mean that funding, competing interests, target-journal formatting,
professional plagiarism/self-plagiarism screening, or DOI-backed raw-data
archiving have been completed.

The 2026-06-20 manuscript pass additionally added the Monte Carlo validation
and uncertainty protocol, the LTE residual mesh diagnostic including the
requested `n=15` case, a timestep sensitivity table, and finite-window current
diagnostics.  The latest path audit checks 34/34 manuscript paths and figures.

## 1. Decisions needed from the authors

Please fill or confirm the following before formal submission.

### Target journal

See `target_journal_shortlist_2026-06-19.md` for a preliminary shortlist based
on official journal/publisher pages checked on 2026-06-19 and spot-verified for
the leading candidates on 2026-06-20.  The current practical recommendation is
SIADS as the first target, with Physica D as the strongest nonlinear-physics
alternative and Journal of Statistical Physics as the strongest
statistical-mechanics alternative.  If the authors choose SIADS, use
`siads_first_submission_packet_2026-06-20.md` for the SIADS-specific cover
letter draft, supplementary-material index, keywords/MSC candidates, and
conversion checklist.

- Target journal:
- Article type:
- Does the journal require a specific LaTeX class/template? yes / no
- Does the journal require anonymized review? yes / no
- Does the journal require ORCID identifiers? yes / no
- Does the journal require data/code DOI rather than GitHub URL? yes / no
- Does the journal require raw data, or are processed/reproducibility artifacts
  sufficient? raw data required / processed artifacts sufficient / unclear

### Author metadata

- Final author order:
- Jayleen Jiang affiliation:
- Jayleen Jiang email:
- Jayleen Jiang ORCID, if any:
- Yao Li affiliation:
- Yao Li email:
- Yao Li ORCID, if any:
- Corresponding author:
- Corresponding-author email:

### Declarations

- Funding:
  - no external funding / funded by the following grants:
  - grant names and numbers:
- Competing interests:
  - none / disclose:
- Ethics:
  - confirm that the current statement is accurate: numerical/theoretical study
    only; no human participants, animal subjects, or personal data.
- AI-assisted preparation:
  - confirm the current disclosure wording or provide journal-specific wording.

### Authorship and contributions

Current manuscript wording:

> Jayleen Jiang performed the numerical experiments, assembled the
> computational artifacts, and drafted the manuscript. Yao Li supervised the
> project and contributed to the model formulation, theoretical framing, and
> interpretation. Both authors should review and approve the final submitted
> version.

Please confirm, edit, or replace this statement.

If the target journal uses CRediT taxonomy, a possible starting point is:

| Contributor | Draft CRediT roles to confirm |
|---|---|
| Jayleen Jiang | Conceptualization; Software; Formal analysis; Investigation; Data curation; Visualization; Writing — original draft; Writing — review and editing |
| Yao Li | Conceptualization; Supervision; Methodology; Writing — review and editing |

## 2. Data/code release decision

Current manuscript state:

- GitHub release bundle is prepared and locally checked.
- `submission_bundle_manifest.md` reports no missing or untracked required
  release files.
- `raw_data_archive_manifest.md` identifies a compact raw-data subset:
  40 unique source-trace raw files, all present locally, totaling
  138,875,181 bytes.
- The full local raw roots are much larger (`Energy Cascade/`, `KDE/`, and
  `lte/`; `lte/` alone is multi-GB), so they should not be committed to GitHub
  without a deliberate release policy.

Choose one:

1. **GitHub-only derived-artifact release.**
   Use the current GitHub branch/release plus the source-trace JSON, claim
   audit, and generated figures.  Suitable only if the journal accepts derived
   reproducibility artifacts and does not require raw training/histogram files.

2. **GitHub + DOI-backed minimal raw-data archive.**
   Archive the 40 files in `raw_data_archive_manifest.md`, preserving their
   `raw_data/...` paths, then replace or supplement the data availability
   statement with the archive DOI.  This is the recommended route if the
   journal asks for raw data.

3. **Full local raw-root archive.**
   Archive the full `Energy Cascade/`, `KDE/`, and `lte/` roots.  This is much
   larger and should be used only if the journal explicitly requires all raw
   historical artifacts rather than the source-trace subset.

Possible DOI-backed data availability wording after an archive exists:

> The manuscript source, analysis scripts, generated figures, and audit
> manifests are available at the project GitHub repository
> [repository URL/release tag].  The raw data files needed to reproduce the
> source-traced numerical tables and figures are archived at [Zenodo/OSF DOI],
> with file paths and SHA-256 checksums listed in
> `raw_data_archive_manifest.md`.

## 3. Professional originality/self-plagiarism check

Local web spot-check status:

- `originality_spotcheck_2026-06-19.md` found no visible exact external phrase
  reuse in sampled web-search snippets.
- This is not a professional similarity report.

Required author action:

- Run iThenticate/Turnitin or the journal/institution-required equivalent on
  the final manuscript.
- Compare against prior drafts, reports, preprints, theses, and any submitted
  manuscripts by the author team.
- Review any overlap with the cited `HLNS` work and any internal project
  reports.

## 4. Generic cover letter template

Replace bracketed text before use.

```text
Dear [Editor Name / Editors],

We are pleased to submit our manuscript entitled
"Numerical study of an energy cascade model derived from a dispersive
equation" for consideration as a [Article Type] in [Journal Name].

The manuscript studies a stochastic energy-cascade chain derived from the
resonant toy model for the cubic nonlinear Schrödinger equation.  The work is a
numerical companion to the rigorous three-mode nonequilibrium steady-state
program and focuses on long-chain behavior that remains analytically open.

The main contributions are:

1. formulation and validation of a Gibbs-preserving two-bath numerical model;
2. high-throughput finite-chain simulations showing an action-current scaling
   E[J(n)] ≈ 28.75 n^{-1.850} over n = 10, 20, 30, 40 under the reported
   parameters;
3. local-equilibrium diagnostics for long-chain pair marginals, including
   explicit limitations where strict local Gibbs structure fails;
4. source-traced short-chain neural Fokker--Planck diagnostics, including
   equilibrium validation, qualitative symmetry breaking, phase locking, and
   slow-mode diagnostics;
5. a reproducibility package with scripted numerical claim audits, path audits,
   release-bundle manifests, and raw-data archive checks.

We believe the manuscript will be of interest to readers working on
nonequilibrium statistical mechanics, stochastic dynamics, wave turbulence, and
numerical studies of Hamiltonian energy cascades.

The manuscript is not under consideration elsewhere.  All authors have approved
the submitted version.  [Add funding, competing-interest, data/code, and
ethical-declaration statements required by the target journal.]

Sincerely,

[Corresponding Author Name]
[Affiliation]
[Email]
```

## 5. Suggested email to coauthor(s)

```text
Subject: Final author/journal confirmations for NLS numerical paper

Hi [Name],

The revised NLS numerical manuscript now has local reproducibility gates in
place.  The command

python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex

currently returns PASS_WITH_LOCAL_RAW_DATA_LIMITATION: the manuscript compiles,
the local numerical claim audit passes 14/14 checks, the path audit passes
34/34 checks, and the release bundle has no missing or untracked required
files.  The remaining limitation is that a DOI-backed raw-data archive has not
yet been uploaded.

Could you please confirm:

1. target journal and article type;
2. final author order, affiliations, corresponding author, and ORCID IDs;
3. funding statement;
4. competing-interest statement;
5. whether the current author-contribution statement is accurate;
6. whether the journal requires raw data in a Zenodo/OSF archive;
7. whether you can run or request the professional similarity check
   (iThenticate/Turnitin/journal equivalent).

Once these items are confirmed, the declarations and data-availability language
can be finalized and the local checks rerun.

Best,
[Your Name]
```

## 6. Final local command after author edits

After inserting author-supplied declarations and any target-journal formatting
changes, rerun:

```sh
python3 Paper/revision_2026-06-19/scripts/run_submission_checks.py --compile-latex
```

Then update or regenerate a final dated integrity-audit snapshot before formal
submission.
