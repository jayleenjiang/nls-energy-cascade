"""Reduce the inputs of Figs. 17-19 to the small arrays in ../results/figure_data/.

Inputs (large, not in the repository):
  MODES   folder with neq_dt1e-3_modes.npz, eq_dt1e-3_modes.npz and the metadata neq_dt1e-3.json, eq_dt1e-3.json
          (written by analyze_step2.py from the sec4_3_spectrum.cpp output, jobs_as_run/run_step2.sbatch)
  SAMPLES folder with n3_neq_val1.bin, n3_eq_val1.bin (independent n = 3 trajectories of Section 4.1,
          ../../sec4_1_stationary_density/jobs_as_run/prep_labels.sbatch)
  EX      folder with the ten files ex_c{0,2}_neq_g{0.1,0.05,0.02,0.01,0.001}.json written by sec4_3_return.cpp
          (jobs_as_run/run_ex.sh, local_runs_g005_g001.sh, runs_g002.sh)
Usage: python3 reduce_figure_data.py MODES SAMPLES EX OUTDIR
Writes
  eigenfunction_coefficients.npz  eigenvalues and the coefficients of psi_0..psi_3 in the dictionary D6 (after fixing
                                  phase, scale and sign on the samples), the dictionary indices and scales
  psi_MH.npz                      Fig. 17: histogram counts and sums of psi_j on the (M, H) grid
  psi1_theta.npz                  Fig. 18: psi_1 on the (theta_1, theta_3) grid at I_1 = I_2 = I_3 = I
  return_probabilities.npz        Fig. 19: mean and standard error of every observable of the ex runs
"""
import json
import os
import struct
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import eigfun as E          # dictionary evaluation and normalisation of the eigenfunctions (uses gedmd_n3.py)

MODES, SAMPLES, EX, OUT = sys.argv[1:5]
os.makedirs(OUT, exist_ok=False)
E.R = MODES.rstrip("/") + "/"


def read_states(path):
    """density-code dump: magic(8) n(u64) d(u32) b(u32), then float32 rows -> (n, d) float64"""
    with open(path, "rb") as f:
        magic = f.read(8); n, d, b = struct.unpack("<QII", f.read(16)); x = np.fromfile(f, "<f4")
    assert magic.startswith(b"NLS") and b == 4 and x.size == n * d, path
    return x.reshape(n, d).astype(np.float64)


# ---- eigenfunctions on 2e5 stationary samples ----------------------------------------------------------------
rng = np.random.default_rng(1)
cases = [("neq_dt1e-3", "n3_neq_val1.bin"), ("eq_dt1e-3", "n3_eq_val1.bin")]
coef, mh, th = {}, {}, {}
N = 121; grid = np.linspace(-np.pi, np.pi, N); T1g, T3g = np.meshgrid(grid, grid)
for case, f in cases:
    lam, C, terms, s = E.load_modes(case)
    xs = read_states(os.path.join(SAMPLES, f)); xs = xs[rng.choice(len(xs), 200000, replace=False)]
    x = np.column_stack([xs[:, :3], -xs[:, 3], xs[:, 4]])          # (I1, I2, I3, theta1, theta3), theta1 = -theta2
    Cj = C[:, :4]; V = E.psi_eval(x, terms, s, Cj)
    Cr, Vr = E.realify(Cj, V, x[:, 1])                              # phase, unit variance, corr(psi, I2) >= 0
    Z = np.load(os.path.join(MODES, case + "_modes.npz"), allow_pickle=True)
    coef[f"{case}_lam"] = lam; coef[f"{case}_C"] = Cr; coef[f"{case}_idx"] = Z["D6_idx"]; coef[f"{case}_scales"] = s
    print(case, np.round(lam[:4], 3))
    # Fig. 17: conditional mean of psi_j given (M, H) = histogram sums / counts on a 40 x 40 grid
    I1, I2, I3, t1, t3 = x.T; M = I1 + I2 + I3
    H = M**2 - 0.5*(I1**2+I2**2+I3**2) + 2*I1*I2*np.cos(t1) + 2*I2*I3*np.cos(t3)
    mb = np.linspace(np.percentile(M, .5), np.percentile(M, 99.5), 41); hb = np.linspace(np.percentile(H, .5), np.percentile(H, 99.5), 41)
    mh[f"{case}_Medges"], mh[f"{case}_Hedges"] = mb, hb
    mh[f"{case}_count"] = np.histogram2d(M, H, [mb, hb])[0]
    for j in (1, 2, 3):
        mh[f"{case}_sum_psi{j}"] = np.histogram2d(M, H, [mb, hb], weights=np.real(Vr[:, j]))[0]
    # Fig. 18: psi_1 at equal actions
    for I in (1.0, 1.5, 2.0):
        xg = np.column_stack([np.full(N * N, I)] * 3 + [T1g.ravel(), T3g.ravel()])
        th[f"{case}_I{I:g}"] = np.real(E.psi_eval(xg, terms, s, Cr[:, [1]])[:, 0]).reshape(N, N)
np.savez_compressed(os.path.join(OUT, "eigenfunction_coefficients.npz"), **coef)
np.savez_compressed(os.path.join(OUT, "psi_MH.npz"), **mh)
np.savez_compressed(os.path.join(OUT, "psi1_theta.npz"), grid=grid, **th)

# ---- Fig. 19: E_x[g(X_t)] averaged over the 32 stream groups ---------------------------------------------------
out = {}; meta = {}
for c in ("c0", "c2"):
    for g in ("0.1", "0.05", "0.02", "0.01", "0.001"):
        name = f"ex_{c}_neq_g{g}"
        d = json.load(open(os.path.join(EX, name + ".json"))); Ex = np.array(d["EX"]); NL = d["nlag"]
        Ex = Ex.reshape(len(Ex), NL, -1) / (d["batches"] * d["lanes"] / d["ngroups"])    # per-group averages
        out[f"{name}_t"] = np.arange(NL) * d["snapshot"]
        out[f"{name}_mean"] = Ex.mean(0); out[f"{name}_se"] = Ex.std(0) / np.sqrt(len(Ex))
        meta[name] = {k: v for k, v in d.items() if k != "EX"}
np.savez_compressed(os.path.join(OUT, "return_probabilities.npz"), **out)
json.dump(meta, open(os.path.join(OUT, "return_probabilities_runs.json"), "w"), indent=1)
print("written to", OUT)
