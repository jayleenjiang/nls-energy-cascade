# Provenance

## Git history

- simulator implementation source state: `3d89659432fac0e512a1cc86fea2b63f8f849762`
- frozen KDE protocol and gates: `0b53e573afdd324fa44b4c7cb50c8bc579d40838`
- low-disk operational change (no scientific parameter changed):
  `094361131337eb19216519c92096ed052a0697c0`
- no-extrapolation unsupported-point handling:
  `0faec919ad58ca8504a6d08a660f1c3492feadfc`
- heterogeneous CSV serialization fix:
  `4e423a965ae076039f7aa3b37fac8b3842b93b1f`

## Source and executable hashes

- `flux/NLS_entropy_ft.cpp`:
  `9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2`
- simulator binary:
  `93c4aa6d046f3c35cd0ea1136091fc44ef1b64baf6237e4663ca9e86b995a156`
- final `analyze_kde_ft.py`:
  `59ddfdb9d530116be169a9ce2651b44580f9c7cc666c7a4f4435e9728dc80e7c`

## Input archives

| case | seed | compressed SHA-256 | decompressed SHA-256 |
|---|---:|---|---|
| driven `(10,2)` | 2026083133 | `5979d35dd335574bc01add81d4c8d97759151a62f4fd6acadc111eb6bdccd7a2` | `4f728b3d0e007d704d90734b0888c00ec05b60f09385c6cfd079f3417d7a088f` |
| equilibrium `T=6` | 2026090106 | `807a4971fbc7ffac1ebad6ad3f8e0a1db83945c18fa38378773f37412954ee22` | `dad1506fad5edda62620905a233f00f01f5522c90b55dd22724d861850b2f4eb` |
| equilibrium `T=10` | 2026090110 | `b5b21983809eb696a77f9de48fa015327085385867a649bf3ce6f276930e81a8` | `83b5658c6b0a7639ed2ce57d4b6ff1761fe603babdda2dbd3b7450760128d1cc` |

The driven decompressed hash exactly matches the accepted original direct
sample recorded before restoration.

## Commands

Exact driven restoration command:

```sh
/Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq \
  sample_n3 10 2 3 8 500 20 7813 0.0005 2026083133 8 \
  /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_total_kde_2026-09-02/raw/driven 1
```

Final analysis command:

```sh
python3 /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_total_kde_2026-09-02/analyze_kde_ft.py \
  --driven /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_total_kde_2026-09-02/raw/driven_blocks.csv.zst \
  --equilibrium-t6 /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/raw/T6_blocks.csv.zst \
  --equilibrium-t10 /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/raw/T10_blocks.csv.zst \
  --output /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_total_kde_2026-09-02/analysis
```

Analysis bootstrap seed: `2026090291`; replicates: `2000` whole-stream
resamples.

## Final report artifacts

- PDF SHA-256:
  `9aeb68f1fd06d56a041508eb83ba037bd04bb6fc61d6b63e1e9cfb68fb7b0684`
- TeX SHA-256:
  `619edef4b8122bd5d0c3b45af19bfe5f7fdbdd0736dad0be21a5c6537bf862ad`
- KDE-accuracy figure SHA-256:
  `e3b5b61a3cc60c0e5d17d2aedb2fe597bfb74185915dc515bff3f9b16f00be50`
- medium-support figure SHA-256:
  `d791989211dc5a03b10c522c5333c89bb7f8728c00961edbcadce4bddda404c9`

The PDF was compiled with TeX Live 2023 via `latexmk`, rendered to five PNG
pages with Poppler, and visually inspected.  The final log has no overfull,
underfull, undefined-reference, or rerun warnings.
