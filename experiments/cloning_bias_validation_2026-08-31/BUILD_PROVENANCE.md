# Build provenance

```text
Apple clang version 17.0.0 (clang-1700.3.19.1)
Target: arm64-apple-darwin24.6.0
Thread model: posix
InstalledDir: /Library/Developer/CommandLineTools/usr/bin
```

- Git HEAD: `76530f438b0b3393b9ffc0043dea02c67b41ae70`
- source Git blob: `0f1160ec03f480582639c2d23f7043b1f6a39260`
- source SHA-256: `bc3c27bf62a45aa879c4a8fd3e4d70fe7bdf1e3bcb0d06c0318826e61e615b6d`
- executable SHA-256: `e07fac4bbca41168114acb8f2505d6d258df73b8bdff3c3ecbad51a2de465088`
- build command: `clang++ -O3 -mcpu=native -std=c++17 -Wall -Wextra -Wpedantic -Xpreprocessor -fopenmp -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp source_archive/NLS_entropy_cloning.cpp -o bin/entropy_cloning_v2`
- self-test command: `bin/entropy_cloning_v2 selftest`
