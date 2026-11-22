#!/usr/bin/env bash
#SBATCH --job-name=b08_build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=slurm_build_%j.log
# Build rk4a_chain_v2 (same file as batch_07, sha256 171c60e5...) and check:
#  CHECK1 n = 4, T1 = Tn = 5: <M> vs exact Gibbs 5.1246;
#  CHECK2 gamma = 0 (no bath), worst start n = 4, t = 1000: relative drift of M (conserved by the Hamiltonian part).
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p bin check
g++ -O3 -fopenmp -std=c++17 rk4a_chain_v2.cpp -o bin/rk4a_chain_v2
{ g++ --version | head -1; sha256sum rk4a_chain_v2.cpp bin/rk4a_chain_v2 tasks_A.tsv tasks_B.tsv; } > BUILD_PROVENANCE.txt
cd check
../bin/rk4a_chain_v2 --n 4 --T1 5 --Tn 5 --dt 1e-3 --burn 50 --measure 300 --streams 64 --threads 2 --seed 11 --out gibbs_n4_T5 > gibbs.log
../bin/rk4a_chain_v2 --n 4 --T1 2 --Tn 8 --gamma 0 --init 0.01,50,50,0.01 --phases 2.094395102393,0.000000000000,2.094395102393,4.188790204786 --dt 1e-4 --burn 0 --measure 1000 --interval 1 --streams 16 --threads 2 --seed 701 --out nobath_n4 > nobath.log
python3 - <<'PY'
import numpy as np
d = {l.split()[0]: l.split()[1:] for l in open("gibbs_n4_T5_summary.txt") if l.strip()}
M = sum(float(d[f"meanI_{j}"][0]) for j in range(1, 5)); se = sum(float(d[f"meanI_{j}"][1]) ** 2 for j in range(1, 5)) ** 0.5
print(f"CHECK1 <M> = {M:.4f} +- {se:.4f}  (exact 5.1246, z = {(M - 5.1246) / se:.2f})  nonfinite_steps = {d['nonfinite_steps'][0]}")
a = np.fromfile("nobath_n4_dump.f32", dtype=np.float32).reshape(1, -1, 16, 7).astype(float)
M = a[..., :4].sum(-1)
print(f"CHECK2 gamma = 0, t <= 1000: max |M/M(0) - 1| = {np.abs(M / 100.02 - 1).max():.2e}  (cloud test: 4.3e-08)")
PY
