# Density training v5 (2026-10-11)

Motivation: with long independent validation trajectories, every n=4 NESS network (action-angle and
Cartesian labels alike) underestimates the bond-angle asymmetry |<sin theta_j>| by 3-10%; the error sits in
the low-action region (lowest two quintiles of P = I_{j-1} I_j), and bounded weak test functions
P sin(theta)/(c+P), P cos(theta)/(c+P) detect it (E[L psi] off by 3-5 SE, trajectories within 1.2 SE).

v5 trainers (train_kde_nn_v5.py for n=4, train_kde_nn_n3_v5.py for n=3) = previous trainers + two options:
  --weak-set bounded   adds those test functions (c = 0.25, 1) to the weak stationarity loss;
  --val-layout streams check-moment SE across 256 independent streams of the long validation runs.
Defaults reproduce the previous trainers.  Labels, queries, Gibbs references: the v2 density run
(../unified_v2/density).  Same epochs, seeds, weights as the v2 density training.

    jv=$(sbatch --parsable val_long.sbatch)
    jt=$(sbatch --parsable --dependency=afterok:$jv train_v5.sbatch)
    sbatch --dependency=afterok:$jt eval_v5.sbatch      # needs the v2 density training finished too
