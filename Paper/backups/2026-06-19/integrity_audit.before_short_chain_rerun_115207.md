# Integrity audit snapshot — 2026-06-19

Scope: `Paper/revision_2026-06-19/draft.tex` after the BibTeX refactor and
short-chain diagnostic tightening.

Verdict: **PASS WITH BLOCKING NOTES for submission**.  Bibliographic existence,
ghost citations, and the main current-scaling claims pass the checks below.
Before journal submission, the LTE table and several short-chain neural-network
diagnostics should be regenerated into machine-readable source metrics, and the
authors must confirm funding / author-contribution declarations.

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
| Full LTE table (`tab:lte`) slopes/R²/temperature comparisons | LTE section | `source_trace_metrics.json`; script `scripts/export_source_trace_metrics.py`; hist/profile source hashes recorded for each row. | VERIFIED for all table rows after aligning `T_kin` entries to the exported pair-averaged profile values. |
| Equilibrium Fokker--Planck validation errors `1.5%, 2.6%, 6.4%`; short-chain density figures | short-chain Fokker--Planck section | `source_trace_metrics.json`; archived notebook outputs from `KDE/4:15_NN/FKE_5d_NLS.ipynb`; figure SHA256 records for `eq_validation.png`, `neq_density.png`, `symmetry_breaking.png`, and `Q1_slices.png`. | VERIFIED for archived equilibrium-error values and figure provenance. Notebook diagnostics without saved output (angular-width ratio, phase-locking peak table, middle-mode balance) are recorded as unverified and no longer used as quantitative manuscript claims. |
| Eigen relaxation diagnostic `lambda_diag≈-0.934` and window sensitivity | eigenfunction section | `eigen_fit_sensitivity.json`; script `scripts/analyze_eigen_fit_windows.py`; source `KDE/NLS_backward_Y_train.txt` | VERIFIED as observable-dependent diagnostic rate, not a spectral-gap claim. |
| Symmetry breaking | stabilization section | Manuscript now uses qualitative wording; old unsupported `15%` claim removed. | ACCEPTABLE as qualitative diagnostic; optional quantitative metric can be added later if recomputed from saved arrays. |

## AI research failure mode checklist

| Mode | Status | Evidence / note |
|---|---|---|
| 1. Implementation bug passing self-review | CLEAR for canonical flux experiment and LTE table export; QUALIFIED for short-chain NN archive | Flux code has frozen hashes, sanitizer smoke run, deterministic threading check, equal-temperature and hot/cold gates. LTE table values are regenerated from hist/profile files. Short-chain NN values are exported from archived notebook outputs; TensorFlow rerun remains a recommended reproducibility supplement. |
| 2. Hallucinated citation | CLEAR | 8/8 references verified; no dangling/orphan citation keys. |
| 3. Hallucinated experimental result | CLEAR for flux/eigen/action-profile/LTE-table claims; QUALIFIED for archived NN figures | Main scaling, eigen diagnostics, action profiles, and LTE table entries are tied to JSON/manifest outputs. Short-chain NN equilibrium-error and symmetry metrics are tied to archived notebook outputs; unarchived diagnostics were removed or downgraded. |
| 4. Shortcut reliance | CLEAR / not applicable | Numerical SDE simulation; no ML benchmark shortcut claim. Neural FP solver is now described diagnostically and validated against Gibbs where possible. |
| 5. Bug reframed as insight | CLEAR for corrected current-scaling narrative | Unsupported GC/entropy and spectral-gap claims were demoted; the paper labels finite-size and diagnostic limitations explicitly. |
| 6. Methodology fabrication | CLEAR for flux and LTE workflows; QUALIFIED for short-chain NN workflow | Flux methods match manifest. LTE table computation is scripted. Short-chain NN workflow is traced to notebooks, model/data hashes, and figure hashes, but a clean non-notebook TensorFlow rerun script remains advisable before final release. |
| 7. Frame-lock | CLEAR WITH LIMITATION | Current framing is explicitly numerical/finite-size and no longer claims asymptotic theorem or full GC test. Remaining limitations are acknowledged. |

## Remaining required actions before calling the paper journal-ready

1. Confirm funding information and author-contribution / competing-interest
   declarations with the authors.
2. Run a professional plagiarism/self-plagiarism check outside this repository
   before formal submission; this audit only performed source and metadata
   integrity checks plus local claim/data traceability.
3. Recommended before final release: archive a TensorFlow environment or
   non-notebook rerun script for the short-chain neural-network diagnostics.
4. Perform the final 100% citation/data/claim audit after author declarations
   and the optional TensorFlow reproducibility supplement are decided.
