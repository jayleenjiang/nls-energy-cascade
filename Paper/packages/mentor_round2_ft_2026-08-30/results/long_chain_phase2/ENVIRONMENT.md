# Phase-II execution environment

- Frozen parent Git commit: `d37f5336ae43b1a1326b281a387f3a346993a000`
- Host: Apple M4, 10 logical CPUs, arm64 macOS Darwin 24.6.0
- Compiler: Apple Clang 17.0.0 (`clang-1700.3.19.1`)
- OpenMP: Homebrew `libomp` using the paths in `build_and_selftest.sh`
- Python: 3.14.5
- NumPy: 2.4.6
- SciPy: 1.17.1
- Matplotlib: 3.11.0

The sampler uses five OpenMP threads per process.  Complementary tilt members
run concurrently, so a production pair may use up to ten logical CPUs.
