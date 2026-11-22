# 4.4 localized test functions (advisor's question, 2026-10-07)

Question: does a localized test function (indicator of a 5D ball, e.g. |x-(1,1,1,0,0)|<0.5) show oscillation?

- report_4_4_localized.pdf / .tex : meeting report (figures fig_loc_backward.png, fig_loc_stationary.png)
- code/local_test.cpp      : stationary C_g(t) of ball indicators, accumulated on the fly (dynamics copied unchanged from sampler_n3.cpp)
- code/local_test_v2.cpp   : same + "ex" mode, E_x[g(X_t)] from a fixed initial point (used for all ex_* runs)
  build: g++ -O3 -march=native -fopenmp -std=c++17 -I<eigen3> local_test_v2.cpp -o local_test_v2
  stat:  ./local_test stat OUT T1 T3 gamma dt burn measure snapshot nlag batches seed threads
  ex:    ./local_test_v2 ex OUT T1 T3 gamma dt 0 0 snapshot nlag batches seed threads I1 I2 I3 th1 th3
- code/run_ex.sh           : the eight E_x runs (x0 = (1,1,1,0,0), (2,1,1,0,0); gamma 0.1 NESS/eq, 0.01, 0.001)
- code/an_stat.py          : C_g(t) + jackknife errors -> data/stat_summary.json
- code/ham_orbit.py        : deterministic (gamma=0) orbit from each center, return times -> data/ham_orbit.json
- code/paper_figs_loc.py   : the two report figures
- data/stat_neq.json, stat_eq.json : raw accumulators (1024 trajectories x 5000 time units, snapshot 0.01, dt 1e-3)
- data/ex_*.json           : raw E_x accumulators (65536 trajectories, t <= 5, snapshot 0.01)

Centers: c0 (1,1,1,0,0), c1 (1,1,1,pi,pi), c2 (2,1,1,0,0), c3 (1,1,2,0,0), c4 (1,1,1,pi/2,-pi/2), c5 (1,2,1,pi,pi); radii 0.5, 1.
