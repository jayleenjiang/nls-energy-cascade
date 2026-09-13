# Validation report

## Production completion

The original four-process pipeline completed the two `dt=2.5e-4` cases and
then exited with status 1 when both `dt=1.25e-4` processes terminated before
writing block CSVs.  That failed attempt is preserved in `pipeline.exitcode`
and the `*.initial_failed.log` files.  The two missing cases were rerun with
exactly the frozen parameters under `RECOVERY_AMENDMENT.md`; the recovery
finished with status 0 at `2026-09-13T11:56:42Z`.

The user-authorized power-policy override affected scheduling only: the
already-running recovery processes were allowed to continue on battery.  No
model, integrator, timestep, seed, sample count, estimator, or gate changed.

## Immutable implementation

- Remote protocol-freeze commit: `5dfe27c0301acfe947c234b75e980bbe0201ef72`
- Actual local launch commit recorded by production:
  `d019ef6e3e2bd0cd149c14738ca73ee13718ffbf`
- Recovery amendment remote commit:
  `103cd95b1055ff02e36e4321dcf4390f1fc3f454`
- Production source SHA-256:
  `98e7f8f5f915c8ce02bd8aa10722025c09fd739184b981961692869c9356c0d3`
- Production binary SHA-256:
  `4c4880d721733897d200f2601690da873f5b226aaa1013df23df2351c1dfe7d1`
- Model parameters: `n=10`, `gamma=0.1`, burn-in `500`, base block
  duration `5`, 128 streams, 31,252 blocks per stream, two threads, middle
  bond 5.

Exact commands are recorded in `COMMANDS.tsv` for the original launch and
`RECOVERY_COMMANDS.tsv` for the two recovered fine-timestep runs.  Seeds are:

| Case | Temperatures | dt | Seed |
|---|---|---:|---:|
| weak_dt2p5e4 | (6.5,5.5) | 2.5e-4 | 2026091101 |
| weak_dt1p25e4 | (6.5,5.5) | 1.25e-4 | 2026091102 |
| moderate_dt2p5e4 | (8,4) | 2.5e-4 | 2026091103 |
| moderate_dt1p25e4 | (8,4) | 1.25e-4 | 2026091104 |

## Raw-data integrity

The independent streaming audit checked the exact CSV header, all numeric
fields for finiteness, all 16,001,024 data rows, and the expected contiguous
mapping `stream_id=floor(row/31252)`, `block_id=row mod 31252`.

| Case | Rows | Streams | SHA-256 |
|---|---:|---:|---|
| weak_dt2p5e4 | 4,000,256 | 128 | `c2c81dc939a8de7f2abf28656b4e5e3140824125586b2db887ea4061243c3f7d` |
| weak_dt1p25e4 | 4,000,256 | 128 | `1b5251958787876634cc484e5c727a008492c33d3e92c3056ee7be5dd7f13192` |
| moderate_dt2p5e4 | 4,000,256 | 128 | `8a7df0b5ec26aba57508702e2ea2552bf7c86b0baeaf4f46e358ef43c7c11b47` |
| moderate_dt1p25e4 | 4,000,256 | 128 | `315bf4a25af37e6d9154ddb423f54cf9af8cd119be9e62f48420e5c17687bc3a` |

The SHA-256 values of all eight inputs used by the analysis, including the two
`dt=5e-4` driven baselines and two equilibrium baselines, were recomputed and
matched `analysis/input_hashes.csv`.

## Frozen analysis checks

- The matched-bin rule, raw support thresholds, two duration selections,
  `1/tau` model, quadratic diagnostic, and 1,000 whole-stream bootstraps are
  unchanged from `PROTOCOL.md`.
- Output counts are exact: 72 per-duration rows, 8,000 per-duration bootstrap
  rows, 16 extrapolation rows, 16,000 extrapolation bootstrap rows, four
  `dt -> 0` rows, and 4,000 `dt -> 0` bootstrap rows.
- Every full-sample resolved per-duration fit and every finite-time
  extrapolation has 1,000 accepted bootstrap replicates.  Nine unsupported
  moderate-affinity long-duration rows remain explicitly unresolved.
- The numerical verdict is
  `MODERATE_AFFINITY_ASYMPTOTE_UNRESOLVED`; see `FINAL_VERDICT.md` and
  `analysis/VERDICT.json` for the strict claim boundary.

## Reproducibility

The full analysis was rerun from the eight raw inputs using the frozen source
and bootstrap seeds.  Every numerical CSV, JSON verdict, and PNG figure was
byte-identical to the pre-rerun artifact.  The four vector PDF figure hashes
changed only because Matplotlib embeds a fresh creation timestamp; their PNG
renders remained byte-identical.  The independently compiled report contains
eight pages and all four figures.

Machine-readable audit details are in `VALIDATION_AUDIT.json`.
