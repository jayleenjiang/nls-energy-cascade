#!/usr/bin/env bash
# Sections 4.1 and 4.2 with the single Cartesian integrator (nls_cartesian v2, co-rotating frame): data step.
#   Gibbs T=5 Metropolis references (n = 3, 4), independent trajectories (1e6 samples each),
#   query states, streamed KDE labels (1e9 samples, h = 0.2, rcut 4, reflection) for
#   n = 3 and n = 4, equilibrium (5,5) and NESS (2,8).  dt = 1e-3, burn 20, interval 0.005.
#SBATCH --job-name=v2prep
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=24G
#SBATCH --time=08:00:00
#SBATCH --output=slurm_prep_%j.log
set -euo pipefail
source /home/jiang37j_mtholyoke_edu/NLS/nls7d_stream_2026-10-03/venv/bin/activate
NCPU=${SLURM_CPUS_PER_TASK:-16}
echo "prep on $(hostname), cpus=$NCPU, $(date)"
mkdir -p eigen && tar xzf eigen-3.4.0-headers.tar.gz -C eigen
g++ -O3 -march=native -fopenmp -std=c++17 -I"$PWD/eigen/eigen3" nls_cartesian.cpp -o nls_cartesian
./nls_cartesian --selftest

python3 gibbs_reference_n.py 2026101000 gibbs_n3_T5.bin 3 > gibbs_n3.log 2>&1
python3 gibbs_reference_n.py 2026101001 gibbs_n4_T5.bin 4 > gibbs_n4.log 2>&1
grep -A4 mean_actions gibbs_n3_T5_diagnostics.json gibbs_n4_T5_diagnostics.json

# independent trajectories: 32 streams x 31250 snapshots = 1e6 samples each
TRAJ="--dt 1e-3 --burn 20 --interval 0.005 --measure 156.25 --batches 2 --threads 2 --checkpoints 1 --dump 1"
for n in 3 4; do
  ./nls_cartesian --n $n --T1 2 --Tn 8 --bc 1 --seed 20261010${n}1 --out n${n}_neq_src $TRAJ > dumps_n$n.log 2>&1 &
  ./nls_cartesian --n $n --T1 2 --Tn 8 --bc 1 --seed 20261010${n}2 --out n${n}_neq_val1 $TRAJ >> dumps_n$n.log 2>&1 &
  ./nls_cartesian --n $n --T1 2 --Tn 8 --bc 1 --seed 20261010${n}3 --out n${n}_neq_val2 $TRAJ >> dumps_n$n.log 2>&1 &
  ./nls_cartesian --n $n --T1 5 --Tn 5 --bc 1 --seed 20261010${n}4 --out n${n}_eq_val1 $TRAJ >> dumps_n$n.log 2>&1 &
  ./nls_cartesian --n $n --T1 5 --Tn 5 --bc 1 --seed 20261010${n}5 --out n${n}_eq_val2 $TRAJ >> dumps_n$n.log 2>&1 &
done
wait
grep -H "mean_I\|mean_current" n?_*_summary.txt

python3 make_queries_n3.py n3_neq_src.bin n3_queries_neq.bin 20000 1
python3 make_queries_n3.py gibbs_n3_T5.bin n3_queries_eq.bin 2000 2
python3 make_queries.py n4_neq_src.bin n4_queries_neq.bin 20000 1
python3 make_queries.py gibbs_n4_T5.bin n4_queries_eq.bin 2000 2

# streamed KDE labels: 256 streams x 3906250 snapshots = 1e9 samples
KDE="--dt 1e-3 --burn 20 --interval 0.005 --measure 19531.25 --batches 16 --threads $NCPU --checkpoints 20 --hI 0.2 --hT 0.2 --rcut 4 --reflect 1"
./nls_cartesian --n 3 --T1 5 --Tn 5 --bc 1 --seed 11 --queries n3_queries_eq.bin  --out n3_eq_kde  $KDE > n3_eq_kde.log 2>&1
python3 eval_eq_n3.py n3_queries_eq.bin n3_eq_kde_kde.bin 5
./nls_cartesian --n 3 --T1 2 --Tn 8 --bc 1 --seed 12 --queries n3_queries_neq.bin --out n3_neq_kde $KDE > n3_neq_kde.log 2>&1
./nls_cartesian --n 4 --T1 5 --Tn 5 --bc 1 --seed 13 --queries n4_queries_eq.bin  --out n4_eq_kde  $KDE > n4_eq_kde.log 2>&1
python3 eval_eq.py n4_queries_eq.bin n4_eq_kde_kde.bin 5
./nls_cartesian --n 4 --T1 2 --Tn 8 --bc 1 --seed 14 --queries n4_queries_neq.bin --out n4_neq_kde $KDE > n4_neq_kde.log 2>&1
grep -H "elapsed\|mean_I\|mean_current\|midpoint_failures" n?_*_kde_summary.txt
echo "prep done $(date)"
