#!/usr/bin/env bash
# Small-gamma worst-case test (2026-10-29), run in the cloud container (2 cores), not on Unity.
# Integrator: rk4a_chain_v2 (same file as batch_07, sha256 171c60e5...). Build: g++ -O3 -fopenmp -std=c++17 rk4a_chain_v2.cpp -o v2
# n = 3..6, gamma = 1e-3, 1e-4, 1e-5, 1e-6; worst (all angles 4pi/3, phases from make_tasks.py) and random phases.
# I_1 = I_n = 0.01, middle modes share M0 = 100 equally; (T1, Tn) = (2, 8); dt 1e-4; t in [0, 0.3]; sampled every 1e-4; 64 trajectories; seed 701.
set -e
for n in 3 4 5 6; do
  mid=$(python3 -c "print(','.join(['%.10g' % (100/($n-2))]*($n-2)))"); init="0.01,$mid,0.01"
  ph=$(python3 -c "
import numpy as np
exec(open('make_tasks.py').read().split('RUN =')[0])
print(phases($n, 4*np.pi/3))")
  for g in 1e-3 1e-4 1e-5 1e-6; do
    A="--n $n --T1 2 --Tn 8 --gamma $g --init $init --dt 1e-4 --burn 0 --measure 0.3 --interval 1e-4 --streams 64 --threads 2 --seed 701"
    ./v2 $A --phases $ph --out n${n}_g${g}_worst
    ./v2 $A --out n${n}_g${g}_random
  done
done
