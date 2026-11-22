#!/usr/bin/env bash
# Build check for batch_09: compile nls_cartesian_v3 and v2 (reference), run the self test, and verify that
#   CHECK1: v3 reproduces the dense BC1 (4,6) n = 20 run of Sec. 3.1 (reference/dense_BC1_T4_T6_n20_summary.txt):
#           mean_I and mean_current must be identical;
#   CHECK2: v2 and v3 give identical summaries on a short run (v3 only adds the line se_current).
#SBATCH --job-name=b09_build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --output=slurm_build_%j.log
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p check/eigen
tar xzf eigen-3.4.0-headers.tar.gz -C check/eigen
g++ -O3 -march=native -fopenmp -std=c++17 -Icheck/eigen/eigen3 nls_cartesian_v3.cpp -o check/v3
g++ -O3 -march=native -fopenmp -std=c++17 -Icheck/eigen/eigen3 nls_cartesian_v2_reference.cpp -o check/v2
{ g++ --version | head -1; sha256sum nls_cartesian_v3.cpp nls_cartesian_v2_reference.cpp tasks_A.csv tasks_B.csv; } > BUILD_PROVENANCE.txt
cd check
./v3 --selftest
./v3 --n 20 --T1 4 --Tn 6 --bc 1 --closed 0 --gamma 0.1 --dt 1e-3 --burn 1280 --measure 2000 --interval 0.1 --batches 4 --threads 4 --checkpoints 4 --seed 202610110084 --out repro_n20 > repro.log
A="--n 12 --T1 4 --Tn 6 --bc 3b --closed 0 --gamma 0.1 --dt 1e-3 --burn 50 --measure 200 --interval 0.1 --batches 3 --threads 3 --checkpoints 4 --seed 77"
./v2 $A --out short_v2 > short_v2.log; ./v3 $A --out short_v3 > short_v3.log
python3 - <<'PY'
rd = lambda f: {l.split()[0]: l.split()[1:] for l in open(f) if l.strip()}
a, b = rd("repro_n20_summary.txt"), rd("../reference/dense_BC1_T4_T6_n20_summary.txt")
ok = all(a[k] == b[k] for k in ("mean_I", "mean_current"))
print("CHECK1 v3 reproduces dense_BC1_T4_T6_n20 (mean_I, mean_current identical):", "yes" if ok else "NO")
s2 = [l for l in open("short_v2_summary.txt") if not l.startswith("elapsed")]
s3 = [l for l in open("short_v3_summary.txt") if not l.startswith(("elapsed", "se_current"))]
print("CHECK2 v2 and v3 summaries identical apart from se_current:", "yes" if s2 == s3 else "NO")
PY
rm -rf eigen
