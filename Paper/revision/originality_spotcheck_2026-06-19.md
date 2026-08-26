# Originality spot-check — 2026-06-19

Scope: local web exact-phrase spot-check for
`Paper/revision_2026-06-19/draft.tex`.  The check used only short manuscript
fragments and author/title queries; the full manuscript was not uploaded to an
external plagiarism service.

Verdict: **PASS WITH LIMITATION**.  The returned web results were topical or
broadly related, but no exact external reuse of the sampled manuscript phrases
was visible in the returned titles/snippets.  Expected public prior work on the
three-mode cascade, especially `HLNS`, did appear and is already cited in the
manuscript.  This spot-check is not a substitute for a professional
plagiarism/self-plagiarism report such as iThenticate or Turnitin before formal
submission.

## Method

- Queried twelve short fragments sampled from the abstract, introduction,
  long-chain method/results, local-equilibrium diagnostics, short-chain
  neural-Fokker--Planck discussion, and conclusion.
- Queried four author/title combinations to detect obvious overlap with public
  preprints, mirrors, or prior project pages.
- Used exact-phrase quotation where practical and inspected the visible result
  titles/snippets for exact external phrase reuse.

## Query log and interpretation

| Group | Query / sampled fragment | Representative returned sources | Interpretation |
|---|---|---|---|
| Abstract / framing | `Gibbs-preserving fluctuation--dissipation bath` | Gibbs-preserving operations pages such as <https://arxiv.org/abs/2404.03479> and <https://link.aps.org/doi/10.1103/PhysRevLett.134.170201> | Topical overlap only; no visible exact phrase reuse. |
| Abstract / scaling claim | `finite-chain action conductivity scales as` | General finite-size thermal-conductivity results, e.g. <https://www.researchgate.net/publication/337180461_Finite-size_effect_of_the_thermal_conductivity_in_one_dimensional_chain> | Topical overlap only; no visible exact phrase reuse. |
| Introduction / model positioning | `forced and dissipated reduction of the nonlinear Schrödinger equation` | Forced NLS literature, e.g. <https://link.aps.org/doi/10.1103/PhysRevE.85.046607> | Related subject matter only; no visible exact phrase reuse. |
| LTE diagnostics | `strict local Gibbs structure fails even at equilibrium` | Broad statistical-mechanics and fluctuation--dissipation material, e.g. <https://en.wikipedia.org/wiki/Fluctuation-dissipation_theorem> | No visible exact phrase reuse. |
| Short-chain validation | `low-accuracy Monte Carlo estimate of the density` | General density-estimation/noise pages, e.g. <https://pmc.ncbi.nlm.nih.gov/articles/PMC12799239/> | No visible exact phrase reuse. |
| Short-chain forcing/noise | `spatially uncorrelated noise in v` | Noise-correlation pages, e.g. <https://www.nature.com/articles/s41534-024-00842-9> and <https://stackoverflow.com/questions/63816481/faster-method-for-creating-spatially-correlated-noise> | Related terminology only; no visible exact phrase reuse. |
| Numerical limitations | `timestep pilot shows that the coarse` | Broad numerical-analysis/noise results | No visible exact phrase reuse. |
| Flux diagnostics | `finite-time averaged currents descriptively` | Broad simulation/statistics results | No visible exact phrase reuse. |
| Short-chain interpretation | `qualitative evidence of symmetry breaking rather than as a standalone` | Symmetry-breaking pages such as <https://link.aps.org/doi/10.1103/PhysRevResearch.2.023244>, <https://arxiv.org/html/2310.02299v7>, and <https://plato.stanford.edu/archives/fall2016/entries/symmetry-breaking/> | Topical overlap only; no visible exact phrase reuse. |
| Eigen diagnostic | `observable-dependent relaxation rate for the observable` | Broad relaxation-rate and eigenfunction results | No visible exact phrase reuse. |
| Limitations / entropy | `bath-energy entropy production and a genuine Gallavotti` | Gallavotti--Cohen and entropy-production literature | Related concepts only; no visible exact phrase reuse. |
| Conclusion / cascade geometry | `open-chain terminal-energy scaling suggested by the cascade geometry` | Broad cascade and finite-chain transport results | No visible exact phrase reuse. |
| Author/title query | `Jayleen Jiang` + `Yao Li` + `energy cascade` | Public `HLNS` prior work, e.g. <https://arxiv.org/abs/2505.16018> and mirror/review pages | Expected related prior work appears; not evidence of uncited overlap. |
| Author/title query | `Jayleen Jiang` + `Numerical study of an energy cascade` | No visible exact public manuscript match in returned titles/snippets | No obvious public duplicate detected. |
| Author/title query | `Yao Li` + `Non-equilibrium steady state` + `three-mode energy cascade` | `HLNS` arXiv and mirrors, e.g. <https://arxiv.org/abs/2505.16018> | Expected cited prior work. |
| Title query | `Non-equilibrium steady state for a three-mode energy cascade model` | `HLNS` arXiv record and mirrors, e.g. <https://arxiv.org/abs/2505.16018> | Confirms public prior work title and citation target. |

## Residual risks

1. The check is a sampled web search, not a corpus-scale plagiarism scan.
2. It cannot reliably detect overlap with paywalled articles, private drafts,
   uploaded course/project reports, or unpublished author material.
3. It does not settle self-plagiarism questions without final author lists,
   prior submissions, and an external similarity report.

Required next step: run a professional plagiarism/self-plagiarism check after
author declarations and journal-format edits are finalized.
