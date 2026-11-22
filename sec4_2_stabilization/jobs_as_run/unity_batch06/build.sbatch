#!/usr/bin/env bash
#SBATCH --job-name=b06_build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=slurm_build_%j.log
# Build rk4a_chain (schemes4 RK4_A + Cartesian bath) and run a short check: n = 4, T1 = Tn = 5, exact Gibbs <M> = 5.1246.
# Built without -march=native so the binary runs on any compute node.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p bin check
g++ -O3 -fopenmp -std=c++17 rk4a_chain.cpp -o bin/rk4a_chain
{ g++ --version | head -1; sha256sum rk4a_chain.cpp bin/rk4a_chain tasks.tsv; } > BUILD_PROVENANCE.txt
cd check
../bin/rk4a_chain --n 4 --T1 5 --Tn 5 --dt 1e-3 --burn 50 --measure 300 --streams 64 --threads 2 --seed 11 --out gibbs_n4_T5
python3 - <<'PY'
d = {l.split()[0]: l.split()[1:] for l in open("gibbs_n4_T5_summary.txt") if l.strip()}
M = sum(float(d[f"meanI_{j}"][0]) for j in range(1, 5)); se = sum(float(d[f"meanI_{j}"][1]) ** 2 for j in range(1, 5)) ** 0.5
print(f"CHECK <M> = {M:.4f} +- {se:.4f}  (exact 5.1246, z = {(M - 5.1246) / se:.2f})  nonfinite_steps = {d['nonfinite_steps'][0]}")
PY
