# Execution environment

- host CPU: Apple M4 (arm64)
- operating system: macOS 15.6.1 (build 24G90)
- compiler: Apple clang 17.0.0 (`clang-1700.3.19.1`)
- C++ standard/optimization: C++17, `-O3 -mcpu=native`
- parallel runtime: Homebrew LLVM OpenMP (`libomp`)
- Python: 3.14.5
- NumPy: 2.4.6
- SciPy: 1.17.1
- Matplotlib: 3.11.0

The production binary is built from `flux/NLS_entropy_cloning.cpp`.  The
committed source manifest records its exact hash; generated binaries are not
versioned.
