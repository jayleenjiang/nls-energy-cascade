# Target-journal policy refresh — SIADS first target

Checked: 2026-06-22.

Scope: official-policy refresh for the current SIADS-first submission plan for
`Paper/revision_2026-06-19/draft.tex` and
`Paper/revision_2026-06-19/draft_siads_review.tex`.  This is a policy and
handoff audit only; it does not change the scientific claims or constitute a
formal submission.

## Official pages checked

Primary official SIADS/SIAM pages:

- SIADS journal page:
  <https://www.siam.org/publications/siam-journals/siam-journal-on-applied-dynamical-systems/>
- SIADS editorial policy:
  <https://epubs.siam.org/journal/siads/editorial-policy>
- SIADS instructions for authors:
  <https://epubs.siam.org/journal/siads/instructions-for-authors>
- SIAM information for journal authors:
  <https://epubs.siam.org/journal-authors>
- SIAM publications editorial policy on artificial intelligence:
  <https://epubs.siam.org/artificial-intelligence>

## Verified policy implications for this manuscript

| Topic | Official-policy implication | Current local status | Action before upload |
|---|---|---|---|
| Journal fit | SIADS publishes mathematical analysis/modeling of dynamical systems and applications, including computational and experimental dynamical-systems work. | Fit remains strong for the stochastic NLS-cascade model, long-chain Monte Carlo, and low-dimensional Fokker--Planck diagnostics. | Keep SIADS as the recommended first target unless the authors prefer a physics/statistical-mechanics audience. |
| Submission files | SIADS submissions use the SIAM Journal Submission & Tracking System; authors provide a cover letter and manuscript PDF, with non-supplemental figures embedded inline. | `draft_siads_review.tex` compiles locally to a 27-page PDF; figures are embedded inline. A placeholder cover-letter template compiles locally. | Replace bracketed cover-letter fields and author declarations before upload. |
| SIAM macros / line numbering | Authors are highly encouraged to use SIAM multimedia macros. If SIAM macros are not used, the instructions require line numbering for review. | The current SIADS review-preparation source uses the non-SIAM-macro line-number fallback via `lineno`; this is a defensible review-preparation route, not a final typesetting conversion. | If time allows, install/check SIAM multimedia macros and convert; otherwise keep the line-numbered fallback and mention any nonstandard formatting only if the submission system asks. |
| Length / file size | SIADS instructions state a general 40-page / 10 MB manuscript expectation, with exceptions justified in the cover letter. | Current audited PDF is 27 pages and well below this practical envelope. | Recheck after final declarations and any SIAM macro conversion. |
| Abstract / keywords / MSC | SIADS requires a one-paragraph abstract under 250 words, plus keywords and Mathematics Subject Classification codes. | Current SIADS source has keywords and MSC candidates. | Authors should confirm the final keywords and MSC list. |
| Reproducibility | SIADS strongly encourages publicly accessible software/data or supplementary materials sufficient to reproduce computational results, including method parameters and post-processing details. | The project has a one-command local gate, claim/path/reference audits, source-bundle builder, raw-data manifest, and source-traced numerical artifacts. | Decide whether GitHub release alone is acceptable or whether to create a DOI-backed raw-data archive for the 42-file, 151,605,557-byte raw subset. |
| Supplementary materials | SIADS encourages supplementary materials; items intended for inclusion with the journal are reviewed with the manuscript and require an index with description and justification. Code should be submitted as a compressed archive, and a TXT index is recommended inside code archives. | `siads_first_submission_packet_2026-06-20.md` already contains a proposed supplementary-material index. | Use the index only if submitting supplements directly; otherwise cite the external GitHub/Zenodo/OSF archive. |
| Original scholarship / similarity | SIAM flags substantial duplication of others' or one's own work as poor scholarship and runs Crossref Similarity Check. | Local web-fragment spot check is clean within sampled scope; no professional similarity report is available locally. | Run iThenticate/Turnitin or the institution/journal equivalent after final author-format edits. |
| AI disclosure | SIAM AI policy v2.0, effective May 2026, requires responsible AI use to be reported in Acknowledgements or Declarations, with a statement that authors assume responsibility for all content. | The manuscript has an AI-assisted-preparation statement, but final wording should match SIAM's current policy and the authors' actual usage. | Add/confirm SIAM-specific AI disclosure before final upload; do not list AI as an author. |
| ORCID / author metadata | SIAM encourages ORCID integration through the submission system. | Author metadata is still incomplete locally. | Confirm author order, affiliations, corresponding author, email, and optional ORCID IDs. |

## Current conclusion

The SIADS-first plan remains scientifically and procedurally reasonable.  The
current local package is strong enough for author/advisor review, but it is not
ready for formal upload until the author-only items are resolved.

No new official-policy finding requires changing the numerical results or
claimed scientific scope.  The main practical additions from the refresh are:

1. keep the SIADS line-numbered review source as a fallback unless authors
   choose to install and use SIAM multimedia macros before submission;
2. make the AI-use statement explicitly SIAM-compatible;
3. treat the raw-data route as a submission decision, with DOI-backed archival
   release preferred for maximal reproducibility;
4. keep professional similarity/self-plagiarism screening as a non-local
   blocker.

## Remaining blockers

- target journal confirmed by authors;
- final author metadata, funding, competing-interest, contribution, and
  AI-use statements supplied by authors;
- professional similarity/self-plagiarism report completed;
- GitHub release and optional DOI-backed raw-data archive route chosen;
- final SIADS PDF and cover-letter PDF regenerated after those decisions;
- one-command local gate rerun after final edits.
