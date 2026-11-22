#!/usr/bin/env bash
# M rise test (2026-10-29): many trajectories, first 0.1 time units, to see the small rise of M directly.
# n = 3..6; gamma = 1e-3 (16384 trajectories) and 1e-4 (65536); worst (all angles 4pi/3) and random phases.
# I_1 = I_n = 0.01, middle modes share M0 = 100; (T1, Tn) = (2, 8); dt 1e-4; sampled every 1e-3; seed 901.
set -e
V=../code/v2
for n in 3 4 5 6; do
  mid=$(python3 -c "print(','.join(['%.10g' % (100/($n-2))]*($n-2)))"); init="0.01,$mid,0.01"
  ph=$(python3 -c "
import numpy as np
exec(open('../code/make_tasks.py').read().split('RUN =')[0])
print(phases($n, 4*np.pi/3))")
  for g in 1e-3 1e-4; do
    S=$([ $g = 1e-3 ] && echo 16384 || echo 65536)
    A="--n $n --T1 2 --Tn 8 --gamma $g --init $init --dt 1e-4 --burn 0 --measure 0.1 --interval 1e-3 --streams $S --threads 2 --seed 901"
    $V $A --phases $ph --out n${n}_g${g}_worst > n${n}_g${g}_worst.log 2>&1
    $V $A --out n${n}_g${g}_random > n${n}_g${g}_random.log 2>&1
    tail -1 n${n}_g${g}_worst.log; tail -1 n${n}_g${g}_random.log
  done
done
