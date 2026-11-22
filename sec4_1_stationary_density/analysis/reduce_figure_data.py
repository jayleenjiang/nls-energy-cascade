"""Reduce the large inputs of Figs. 9-12 to the small arrays stored in ../results/figure_data/.

Needs TensorFlow (2.21 was used) and three large inputs that are not in the repository:
  GIBBS  folder with gibbs_n3_T5.bin, gibbs_n4_T5.bin   (Metropolis samples of the T = 5 Gibbs measure, made by
                                                          ../labels/gibbs_reference_n.py in ../jobs_as_run/prep_labels.sbatch)
  VAL    folder with n{3,4}_neq_val{A,B}.bin             (long validation trajectories, ../jobs_as_run/val_long.sbatch)
The trained networks are in ../models.

Usage: python3 reduce_figure_data.py GIBBS VAL OUTDIR
Writes OUTDIR/marginals_n3.npz, marginals_n4.npz (2-d histograms) and OUTDIR/theta_slices.npz (density sections).
The computation is the one of an_sec4_v5.py; ../figures/make_figures.py only draws these arrays.
"""
import os
import struct
import sys

import numpy as np
import tensorflow as tf
from scipy.special import logsumexp

tf.keras.backend.set_floatx("float64")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "training"))
import n3_common, n3_model, n4_common, n4_model   # network definition and reduced Hamiltonian, as used in training

GIBBS, VAL, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
MODELS = os.path.join(HERE, "..", "models")
os.makedirs(OUT, exist_ok=False)            # never write into an existing folder


def read_states(path):
    """density-code dump: magic(8) n(u64) d(u32) b(u32), then float32 rows -> (n, d) float64"""
    with open(path, "rb") as f:
        magic = f.read(8); n, d, b = struct.unpack("<QII", f.read(16)); x = np.fromfile(f, "<f4")
    assert magic.startswith(b"NLS") and b == 4 and x.size == n * d, path
    return x.reshape(n, d).astype(np.float64)


def logratio(m, x, bs=65536):
    """network output log(rho / rho_Gibbs) at the states x, in batches"""
    return np.concatenate([m.log_ratio(tf.constant(x[j:j + bs])).numpy() for j in range(0, len(x), bs)])


# ---- Figs. 9-10: two-dimensional marginals ------------------------------------------------------------------
# The network density is represented by reweighting the Gibbs samples with w ~ exp(log_ratio); the simulation
# density by the histograms of the two independent validation sets A and B.
for n in (3, 4):
    M = n3_model if n == 3 else n4_model
    ref = read_states(f"{GIBBS}/gibbs_n{n}_T5.bin")
    m = M.load_bundle(f"{MODELS}/n{n}_neq_s1")
    f = logratio(m, ref); w = np.exp(f - logsumexp(f))
    A = read_states(f"{VAL}/n{n}_neq_valA.bin"); B = read_states(f"{VAL}/n{n}_neq_valB.bin")
    if n == 3:     # stored angle theta_2 = -theta_1 of the paper
        conv = lambda x: np.column_stack([x[:, :3], -x[:, 3], x[:, 4]])
        rows = [(0, 2), (3, 4), (0, 3)]          # (I1,I3), (theta1,theta3), (I1,theta1)
    else:
        conv = lambda x: x
        rows = [(0, 3), (4, 5), (0, 4)]          # (I1,I4), (theta2,theta3), (I1,theta2)
    refp, Ap, Bp = conv(ref), conv(A), conv(B)
    edges = lambda j: np.linspace(0, 4.5, 19) if j < n else np.linspace(-np.pi, np.pi, 19)
    h2 = lambda x, a, b, w=None: np.histogram2d(x[:, a], x[:, b], [edges(a), edges(b)], weights=w)[0]
    h = len(ref) // 2; wa, wb = w[:h] / w[:h].sum(), w[h:] / w[h:].sum()   # two halves: network noise level
    out = {"rows": np.array(rows)}
    for r, (a, b) in enumerate(rows):
        out[f"r{r}_edges_x"], out[f"r{r}_edges_y"] = edges(a), edges(b)
        out[f"r{r}_simA"], out[f"r{r}_simB"] = h2(Ap, a, b), h2(Bp, a, b)            # counts of the two validation sets
        out[f"r{r}_net"] = h2(refp, a, b, w)                                         # network, all Gibbs samples
        out[f"r{r}_netA"], out[f"r{r}_netB"] = h2(refp[:h], a, b, wa), h2(refp[h:], a, b, wb)
    np.savez_compressed(os.path.join(OUT, f"marginals_n{n}.npz"), **out)
    print("written", f"marginals_n{n}.npz")

# ---- Figs. 11-12: sections of the density in the angles at equal actions I_j = I ------------------------------
g = np.linspace(-np.pi, np.pi, 121); dA = (g[1] - g[0]) ** 2
store = {}
for n in (3, 4):
    M = n3_model if n == 3 else n4_model
    E_np = n3_common.energy_np if n == 3 else n4_common.energy_np
    m = M.load_bundle(f"{MODELS}/n{n}_neq_s1")
    for I in (0.3, 0.5, 1.0):
        if n == 3:
            G1, G3 = np.meshgrid(g, g, indexing="ij")
            xp = np.column_stack([np.full(G1.size, I)] * 3 + [G1.ravel(), G3.ravel()]); xf = xp.copy(); xf[:, 3] = -xf[:, 3]
            E = E_np(xf); fx = logratio(m, xf)
            dens = [np.exp(Lg - Lg.max()) for Lg in ((-E / 5 + fx).reshape(G1.shape), (-E / 5).reshape(G1.shape))]
        else:      # n = 4: average over the third angle theta_4
            t4 = np.linspace(-np.pi, np.pi, 48, endpoint=False)
            G2, G3, G4 = np.meshgrid(g, g, t4, indexing="ij")
            xf = np.column_stack([np.full(G2.size, I)] * 4 + [G2.ravel(), G3.ravel(), G4.ravel()])
            E = E_np(xf); fx = logratio(m, xf); dens = []
            for Lg in ((-E / 5 + fx), (-E / 5)):
                Lg = Lg.reshape(G2.shape); dens.append(np.exp(Lg - Lg.max()).mean(-1))
        for r, p in enumerate(dens):    # r = 0: NESS network density, r = 1: Gibbs density at T = 5
            store[f"n{n}_{'nn' if r == 0 else 'gibbs'}_I{I}"] = p / (p.sum() * dA)
    p = store[f"n{n}_nn_I0.5"]; kk = np.unravel_index(p.argmin(), p.shape)
    print(n, "NN I=0.5 min at", round(g[kk[0]], 2), round(g[kk[1]], 2), "ratio", round(float(p.max() / p.min()), 3))
np.savez_compressed(os.path.join(OUT, "theta_slices.npz"), grid=g, **store)
print("written theta_slices.npz")
