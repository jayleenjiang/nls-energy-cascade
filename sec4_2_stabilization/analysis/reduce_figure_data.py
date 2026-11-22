"""Reduce the trajectory dumps behind Figs. 13-16 to the small arrays in ../results/figure_data/.

Inputs (large, not in the repository; each made by sec4_2_stabilization.cpp with the commands in ../jobs_as_run):
  ANGLES   folder angles/ of jobs_as_run/local_runs/run_rk4a_all.sh   (Fig. 13, middle modes held at 0)
  WORST    worst_case_curves.npz written by worstcase/worst_case_analyze.py on the batch_07 output   (Fig. 14)
  SMALLG   folder of jobs_as_run/local_smallgamma/run_smallg.sh        (Fig. 15)
  MRISE    folder of jobs_as_run/local_Mrise/run_Mrise.sh              (Fig. 16)
Usage: python3 reduce_figure_data.py ANGLES WORST SMALLG MRISE OUTDIR
Dump layout: float32 rows (I_1..I_n, theta_2..theta_n), ordered (stream block of 16, sample, lane).
"""
import os
import sys

import numpy as np

ANG, WORST, SMALLG, MRISE, OUT = sys.argv[1:6]
os.makedirs(OUT, exist_ok=False)

# ---- Fig. 13: stationary densities of the end angles, 72 bins on (-pi, pi) ------------------------------------
out = {}
for case, n in (("n3_hold2", 3), ("n4_hold2", 4), ("n4_hold23", 4)):
    for T in ("2_8", "5_5", "8_2"):
        D = 2 * n - 1
        a = np.fromfile(f"{ANG}/{case}_T{T}_dump.f32", dtype=np.float32).reshape(4, -1, 16, D).reshape(-1, D)
        curves = {"theta1": -a[:, n], "theta3": a[:, n + 1]}      # theta_1 = -theta_2 (end bond oriented outward)
        if n == 4:
            curves["theta4"] = a[:, n + 2]
        for lab, x in curves.items():
            h, e = np.histogram(np.remainder(x + np.pi, 2 * np.pi) - np.pi, bins=72, range=(-np.pi, np.pi), density=True)
            out[f"{case}_T{T}_{lab}"] = h
out["edges"] = e
np.savez_compressed(os.path.join(OUT, "boundary_angles.npz"), **out)

# ---- Fig. 14: mean M(t) for M0 = 200 and the three initial phases -------------------------------------------
Z = np.load(WORST)
np.savez_compressed(os.path.join(OUT, "worstcase_M200.npz"),
                    **{k: Z[k] for k in Z.files if "_M200_" in k and k.endswith(("__t", "__mean"))})

# ---- Fig. 15: median over 64 trajectories of I_1 + I_n --------------------------------------------------------
out = {}
for n in (3, 4, 5, 6):
    d = 2 * n - 1
    for g in ("1e-3", "1e-4", "1e-5", "1e-6"):
        for ph in ("worst", "random"):
            a = np.fromfile(f"{SMALLG}/n{n}_g{g}_{ph}_dump.f32", dtype=np.float32).reshape(4, -1, 16, d).transpose(0, 2, 1, 3).reshape(64, -1, d)
            out[f"n{n}_g{g}_{ph}"] = np.median(a[:, :, 0] + a[:, :, n - 1], 0)    # sampled every 1e-4 from t = 0
np.savez_compressed(os.path.join(OUT, "smallgamma_endpoints.npz"), **out)

# ---- Fig. 16: (E M(t) - M(0)) / gamma with its standard error ------------------------------------------------
out = {}; M0 = 100.02                                   # 0.01 + 100 + 0.01
for g in ("1e-3", "1e-4"):
    for n in (3, 4, 5, 6):
        d = 2 * n - 1
        for ph in ("worst", "random"):
            a = np.fromfile(f"{MRISE}/n{n}_g{g}_{ph}_dump.f32", dtype=np.float32).reshape(-1, 100, 16, d)
            M = a[..., :n].astype(float).sum(-1)
            dM = (M - M0).transpose(1, 0, 2).reshape(M.shape[1], -1)              # (sample time, trajectory)
            out[f"n{n}_g{g}_{ph}_mean"] = dM.mean(1) / float(g)
            out[f"n{n}_g{g}_{ph}_se"] = dM.std(1) / np.sqrt(dM.shape[1]) / float(g)
np.savez_compressed(os.path.join(OUT, "smallgamma_Mrise.npz"), **out)
print("written to", OUT)
