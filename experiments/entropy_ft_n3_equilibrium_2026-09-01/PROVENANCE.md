# Production provenance

- simulator implementation commit: `70ebf16`
- source snapshot used for compilation: branch state `3d89659`
- frozen protocol commit present before output: `c7e2a92`
- operational AC-resume-only commit: `f917173`
- source: `flux/NLS_entropy_ft.cpp`
- source SHA-256: `9ae5835ed708c8794c8b00ba799b23761482953aaf0ed47cd0b4ba3966d4eaf2`
- binary SHA-256: `93c4aa6d046f3c35cd0ea1136091fc44ef1b64baf6237e4663ca9e86b995a156`
- compiler: Apple clang 17.0.0 (`clang-1700.3.19.1`)
- hardware: Apple M4

Build command:

```sh
clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic \
  -Xpreprocessor -fopenmp -I/opt/homebrew/include/eigen3 \
  -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp \
  flux/NLS_entropy_ft.cpp \
  -o experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq
```

Exact production commands:

```sh
/Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq sample_n3 6 6 3 8 500 20 7813 0.0005 2026090106 8 /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/raw/T6 1
/Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/bin/entropy_ft_n3_eq sample_n3 10 10 3 8 500 20 7813 0.0005 2026090110 8 /Users/jayleenjiang/NLS_eq_runtime/experiments/entropy_ft_n3_equilibrium_2026-09-01/raw/T10 1
```

The unchanged sampler wrote each blocks CSV through a FIFO into Zstandard.
Both archives pass `zstd -t`; each contains 1,000,065 decompressed lines.

- `T6_blocks.csv.zst`: 120,749,526 bytes; compressed SHA-256
  `807a4971fbc7ffac1ebad6ad3f8e0a1db83945c18fa38378773f37412954ee22`;
  decompressed CSV SHA-256
  `dad1506fad5edda62620905a233f00f01f5522c90b55dd22724d861850b2f4eb`.
- `T10_blocks.csv.zst`: 120,631,709 bytes; compressed SHA-256
  `b5b21983809eb696a77f9de48fa015327085385867a649bf3ce6f276930e81a8`;
  decompressed CSV SHA-256
  `83b5658c6b0a7639ed2ce57d4b6ff1761fe603babdda2dbd3b7450760128d1cc`.

Analysis command:

```sh
python3 experiments/entropy_ft_n3_equilibrium_2026-09-01/analyze_equilibrium.py \
  experiments/entropy_ft_n3_equilibrium_2026-09-01 \
  experiments/entropy_ft_n3_equilibrium_2026-09-01/analysis
```
