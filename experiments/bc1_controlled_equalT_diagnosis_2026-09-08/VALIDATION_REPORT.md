# Validation report

## Integrity

- 512 trajectories at n=25 and 512 at n=100.
- 200001 measurement samples per trajectory.
- Missing products: none. Integrity errors: none.
- All 426 saved-output SHA-256 entries pass `shasum -a 256 -c`.
- Fresh reanalysis matches `analysis/final_analysis_stdout.json` byte-for-byte.

## Frozen numerical protocol

- BC1, T1=Tn=6, gamma=0.1, dt=5e-4.
- Burn-in: 2000 (n=25), 32000 (n=100); measurement duration 2000; sample interval 0.01.
- Seeds: 2026090825, 2026090826, 2026090900, 2026090901.
- Source SHA-256: `17d2b0e75c7838c7fc8710bae0f0381c7607c967d38e03b4db3bde34e8ec086b`.
- Binary SHA-256: `68c13863da721f655ecbd3d0913b80691c0a0a3a6c984877af8423c331c443db`.
- Launch commit: `908a7c08faac7ea3b8eb06c0577fca7fb5e9ad5e`.
- Remote frozen commit: `35124a825a4119e9007c635e739ae71e3752d6b7`.

## Gates

| n | action port | numerical events | max |z| | pointwise failures | Bonferroni | final control |
|---:|:---:|:---:|---:|---:|:---:|:---:|
| 25 | PASS | PASS | 1.875311 | 0/24 | PASS | PASS |
| 100 | PASS | PASS | 2.984038 | 5/99 | PASS | PASS |

Full per-bond and per-site raw analysis values are in `analysis/*_comparison.csv` and in the PDF appendices. The 64 unique completed commands are in `provenance/COMMANDS_UNIQUE_PRODUCTION.csv`; all 72 launch records, including the eight no-data detached-attempt records, remain in `provenance/COMMANDS.csv`. Frozen output hashes are in `provenance/OUTPUT_HASHES.sha256`.
