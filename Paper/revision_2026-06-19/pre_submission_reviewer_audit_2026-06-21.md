# Pre-submission reviewer-style audit — 2026-06-21

Scope: independent local reviewer-style audit of
`Paper/revision_2026-06-19/draft.tex` after the `n=60` fine-timestep and
journal-positioning updates. This is a manuscript-readiness review, not an
author declaration, plagiarism-screening, or target-journal formatting
certification.

## Reviewer configuration

| Role | Review angle |
|---|---|
| Editor-style reviewer | Journal fit, scope control, and whether the paper tells one coherent story. |
| Methodology reviewer | Monte Carlo design, uncertainty claims, validation gates, and reproducibility. |
| Domain reviewer | Relationship to resonant NLS, wave turbulence, and nonequilibrium transport. |
| Numerical-PDE reviewer | Role and limitations of the neural Fokker--Planck and eigenfunction diagnostics. |
| Devil's advocate | Places where a reader could misread a finite-size or diagnostic statement as a theorem. |

## Editorial decision

Local decision: **minor revision / locally submission-ready after author and
journal confirmations**.

The manuscript now has a defensible finite-size claim structure: the primary
transport exponent remains tied to the matched `n=10,20,30,40` design, while
`n=50`, `n=60`, timestep, bath-temperature, and thermostat-coupling studies are
framed as robustness checks. The claim-evidence map, validation appendix,
source-traced LTE residual diagnostics, and one-command gate make the numerical
scope unusually transparent for a first submission.

## Reviewer consensus

1. **No remaining local numerical blocker was found.** The latest one-command
   gate reports `PASS_WITH_AUTHOR_CONFIRMATION_PENDING_AND_LOCAL_RAW_DATA_LIMITATION`,
   with 22/22 registered numerical claims verified and 11/11 cited BibTeX
   entries structurally checked.
2. **The main scientific vulnerability is correctly disclosed.** The exponent
   is a finite-size law, not an asymptotic theorem. The manuscript explicitly
   says that long-chain ergodicity, asymptotic scaling, bath-energy entropy
   production, and Gallavotti--Cohen large deviations remain open.
3. **The long-chain and short-chain parts now have a coherent division of
   labor.** Long-chain Monte Carlo supports transport/LTE claims; short-chain
   neural Fokker--Planck computations provide mechanistic diagnostics rather
   than proof of the long-chain exponent.
4. **Terminology needed one final presentation fix.** The visible title
   `Thermal conductivity` could invite a reviewer to think the paper claims
   Hamiltonian heat conductivity. This audit therefore recommends using
   `Action-current conductivity` in visible headings while preserving the
   existing LaTeX label for cross-reference stability.

## Action items from this audit

| Priority | Item | Status |
|---|---|---|
| P0 | Keep author declarations, funding, final author order, ORCID, and competing-interest wording pending author confirmation. | External/author item; not edited locally. |
| P0 | Keep professional similarity screening and DOI-backed raw-data upload outside local claims until actually completed. | External item; still pending. |
| P1 | Rename visible conductivity terminology from thermal/heat wording to action-current wording where the measured observable is `J_j`. | Applied in the manuscript sources after this audit. |
| P1 | Keep the old `sec:thermal-conductivity` LaTeX label unless a full reference migration is needed. | Applied; cross-references remain stable. |
| P1 | Register this reviewer-style audit in the release manifest so Gate E has explicit evidence. | Applied via `build_submission_bundle_manifest.py`. |
| P2 | Optional future strengthening: another still-larger chain length, multi-length timestep convergence, or full neural-network retraining. | Not required for the current finite-size claim scope. |

## Residual risk statement

The paper can be handed to the advisor as locally vetted for scientific
consistency and reproducibility. It should not yet be called formally
submitted-ready until the author-supplied declarations, target-journal
formatting decision, professional similarity report, and raw-data archive route
are finalized.
