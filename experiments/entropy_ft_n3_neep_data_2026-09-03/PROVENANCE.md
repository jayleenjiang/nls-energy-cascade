# Production provenance

## Frozen implementation

- validated source snapshot commit:
  `3d89659432fac0e512a1cc86fea2b63f8f849762`
- validated equilibrium-audit completion commit:
  `a4bd3b596ca557ec07b0323012d8bc626089cd30`
- NEEP data protocol commit, present before production:
  `2b945d9b8a38fb24a8c191f6ea168866f80c48bf`
- source: `flux/NLS_entropy_ft.cpp`
- source SHA-256:
  `9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2`
- reused validated binary:
  `experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq`
- binary SHA-256:
  `93c4aa6d046f3c35cd0ea1136091fc44ef1b64baf6237e4663ca9e86b995a156`
- periodic-column converter SHA-256:
  `82c34bca3ae23b2a23b035f431541bd1f0d638f5a587fda385ea6e253e81552c`
- sanity-analysis script SHA-256:
  `45bf222077b331b3c19eb86daeffd50343b65253f501c185e4d66624883f8ce5`

Immediately before production, the unchanged binary self-test reported:

```text
selftest maximum gradient error=1.55486956643e-10
selftest Hamiltonian dE/dt=-9.02131476493e-17
selftest maximum boundary Laplacian error=0.000529314434244
SELFTEST PASS
```

## Exact simulator commands

Driven seed: `2026090310`.

```sh
/Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq \
  sample_n3 10 2 3 8 500 0.1 39063 0.0005 2026090310 8 \
  /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_neep_data_2026-09-03/raw/driven 1
```

Equilibrium seed: `2026090306`.

```sh
/Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq \
  sample_n3 6 6 3 8 500 0.1 39063 0.0005 2026090306 8 \
  /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_neep_data_2026-09-03/raw/equilibrium 1
```

Each unchanged blocks CSV was written through a FIFO and compressed with
`zstd -T0 -12`.  The deterministic NEEP conversion command for each label was
equivalent to:

```sh
zstd -dc raw/LABEL_blocks.csv.zst \
  | python3 augment_neep.py \
  | zstd -T0 -12 -q -c > raw/LABEL_neep_transitions.csv.zst
```

No uncompressed full dataset was materialized.

## Production observations

- driven: 5,000,064 rows, zero midpoint failures, mean midpoint iterations
  `5.39825`, simulator elapsed time `30.331 s`;
- equilibrium: 5,000,064 rows, zero midpoint failures, mean midpoint
  iterations `5.32911`, simulator elapsed time `33.1869 s`.

Both conversion logs report exactly `augmented_rows=5000064`.

## Output hashes

| case | archive | compressed SHA-256 | decompressed CSV SHA-256 |
|---|---|---|---|
| driven | raw | `555ea347d1aa0043d4e10d8aeaa6f22179bdd3e603555b102e8184f2ad46a3c3` | `c8ad530f0475958a22ac062a3059a41b8416b54b8dda2a95a3a53abdb5118fd2` |
| driven | NEEP-ready | `95b7c201eb3da3d1cb3dd0e3a988c9ee043b680fe59510c652a2ecc09eeeda48` | `c9a3edaf93d65c9aa5b4563b876ab3b5f6803bf199da2c638ae2cc613c7af6a4` |
| equilibrium | raw | `b590eb4c56c693bfff37e6bfef30c21517143767d4d38e314eb8eccb60198887` | `b9c1e5c73a9e77256095f87230539a1070e7f6e83f73677ea1f3a0463ce655e7` |
| equilibrium | NEEP-ready | `fb63862ce01809ef2e2d05c0a35797d629b76a2dbe783d33c327a2e8d4a39738` | `0de7c9cfc05fef1b6493b4ef0a586b8b45be5a47725dced08c450b28773d2135` |

The raw production data are intentionally excluded from Git because the four
archives total 2.5655 GiB.  Git stores the frozen protocol, generation and
analysis code, complete hashes, manifest, and numerical audit only.
