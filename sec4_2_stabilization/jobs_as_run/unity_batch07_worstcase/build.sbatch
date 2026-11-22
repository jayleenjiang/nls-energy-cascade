#!/usr/bin/env bash
#SBATCH --job-name=b07_build
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=00:20:00
#SBATCH --output=slurm_build_%j.log
# Build rk4a_chain_v2 (= rk4a_chain of batch_06 + option --phases) and check:
#  (1) without --phases v2 gives byte-identical output to v1 (rk4a_chain_v1_reference.cpp = batch_06 rk4a_chain.cpp);
#  (2) n = 4, T1 = Tn = 5: <M> vs exact Gibbs 5.1246.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
mkdir -p bin check
g++ -O3 -fopenmp -std=c++17 rk4a_chain_v2.cpp -o bin/rk4a_chain_v2
g++ -O3 -fopenmp -std=c++17 rk4a_chain_v1_reference.cpp -o bin/rk4a_chain_v1
{ g++ --version | head -1; sha256sum rk4a_chain_v2.cpp rk4a_chain_v1_reference.cpp bin/rk4a_chain_v2 bin/rk4a_chain_v1 tasks.tsv; } > BUILD_PROVENANCE.txt
cd check
A="--n 5 --T1 2 --Tn 8 --init 0.01,30,30,30,0.01 --dt 1e-4 --burn 0 --measure 2 --interval 0.05 --streams 32 --threads 2 --seed 7"
../bin/rk4a_chain_v1 $A --out v1 > v1.log; ../bin/rk4a_chain_v2 $A --out v2 > v2.log
if cmp -s v1_dump.f32 v2_dump.f32; then echo "CHECK1 v1 and v2 dumps identical: yes"; else echo "CHECK1 v1 and v2 dumps identical: NO"; fi
../bin/rk4a_chain_v2 --n 4 --T1 5 --Tn 5 --dt 1e-3 --burn 50 --measure 300 --streams 64 --threads 2 --seed 11 --out gibbs_n4_T5 > gibbs.log
python3 - <<'PY'
d = {l.split()[0]: l.split()[1:] for l in open("gibbs_n4_T5_summary.txt") if l.strip()}
M = sum(float(d[f"meanI_{j}"][0]) for j in range(1, 5)); se = sum(float(d[f"meanI_{j}"][1]) ** 2 for j in range(1, 5)) ** 0.5
print(f"CHECK2 <M> = {M:.4f} +- {se:.4f}  (exact 5.1246, z = {(M - 5.1246) / se:.2f})  nonfinite_steps = {d['nonfinite_steps'][0]}")
PY
