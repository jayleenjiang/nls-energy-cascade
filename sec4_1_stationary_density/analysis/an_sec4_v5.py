#!/usr/bin/env python3
"""Sections 4.1/4.2 (rev7): numbers, marginal tables and figures from the v5 density models.
Only Cartesian implicit-midpoint data (nls_cartesian v2, co-rotating frame):
  models      : batch_04_density_v5/models/{n3_eq,n3_neq_s1,n3_neq_s2,n4_eq,n4_neq_s1,n4_neq_s2}
  Gibbs ref   : unified_v2/density/gibbs_n{3,4}_T5.bin (Metropolis, 6.4e5 states)
  validation  : long trajectories valA/valB (256 streams, burn 500, measure 1e4 every 0.5, dt 1e-3)
State in files: (I_1..I_n, theta_2..theta_n), theta_j = 2(phi_j - phi_{j-1}).
n = 3 figures use the paper angle theta_1 = -theta_2.
usage: python3 an_sec4_v5.py MODELS_DIR DENSITY_DIR VAL_DIR OUTDIR
"""
import sys, json, struct, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
import tensorflow as tf
tf.keras.backend.set_floatx("float64")
from scipy.special import logsumexp
MOD, DEN, VAL, OUT = sys.argv[1:5]
sys.path.insert(0, ".")
import n3_model, n4_model, n3_common, n4_common

def read_states(path):
    with open(path, "rb") as f:
        magic = f.read(8); n, d, b = struct.unpack("<QII", f.read(16)); x = np.fromfile(f, "<f4")
    assert magic.startswith(b"NLS") and b == 4 and x.size == n * d, path
    return x.reshape(n, d).astype(np.float64)

def logratio(m, x, bs=65536):
    return np.concatenate([m.log_ratio(tf.constant(x[j:j + bs])).numpy() for j in range(0, len(x), bs)])

def currents(x, n):
    I, th = x[:, :n], x[:, n:]
    return np.column_stack([4 * I[:, k] * I[:, k + 1] * np.sin(th[:, k]) for k in range(n - 1)])

res = {}
for n in (3, 4):
    M = n3_model if n == 3 else n4_model
    E_np = n3_common.energy_np if n == 3 else n4_common.energy_np
    ref = read_states(f"{DEN}/gibbs_n{n}_T5.bin")
    obs = np.column_stack([ref[:, :n], currents(ref, n)])
    nb = 50; blk = np.array_split(np.arange(len(ref)), nb)
    R = res[n] = {}
    W = {}
    for case in ("eq", "neq_s1", "neq_s2"):
        m = M.load_bundle(f"{MOD}/n{n}_{case}")
        f = logratio(m, ref); w = np.exp(f - logsumexp(f)); W[case] = w
        mean = w @ obs
        bm = np.array([(w[b] @ obs[b]) / w[b].sum() for b in blk]); se = bm.std(0, ddof=1) / np.sqrt(nb)
        d = dict(profile=mean[:n].tolist(), currents=mean[n:].tolist(), se=se.tolist(),
                 ess=float(1 / np.sum(w * w) / len(w)))
        if case == "eq":
            d["u_centred_rms_gibbs"] = float(np.sqrt(np.mean((f - f.mean()) ** 2)))
        R[case] = d
        print(n, case, {k: (np.round(v, 4).tolist() if isinstance(v, list) else round(v, 4)) for k, v in d.items()}, flush=True)
    # ---- marginals against the long validation trajectories (NESS)
    A = read_states(f"{VAL}/n{n}_neq_valA.bin"); B = read_states(f"{VAL}/n{n}_neq_valB.bin")
    R["val_profile_currents"] = np.concatenate([np.vstack([A, B])[:, :n].mean(0), currents(np.vstack([A, B]), n).mean(0)]).tolist()
    R["val_samples"] = [len(A), len(B)]
    if n == 3:
        conv = lambda x: np.column_stack([x[:, :3], -x[:, 3], x[:, 4]])     # paper (I1,I2,I3,theta1,theta3)
        lab = ["I₁", "I₂", "I₃", "θ₁", "θ₃"]; planes = [(0, 1), (0, 2), (3, 4), (0, 3), (2, 4)]; rows = [(0, 2), (3, 4), (0, 3)]
    else:
        conv = lambda x: x
        lab = ["I₁", "I₂", "I₃", "I₄", "θ₂", "θ₃", "θ₄"]; planes = [(0, 1), (0, 3), (4, 5), (0, 4), (3, 6)]; rows = [(0, 3), (4, 5), (0, 4)]
    refp, Ap, Bp = conv(ref), conv(A), conv(B)
    def edges(j): return np.linspace(0, 4.5, 19) if j < n else np.linspace(-np.pi, np.pi, 19)
    def ext(j): return (0, 4.5) if j < n else (-np.pi, np.pi)
    def h2(x, a, b, w=None): return np.histogram2d(x[:, a], x[:, b], [edges(a), edges(b)], weights=w)[0]
    h = len(ref) // 2; tab = {}; figrows = []
    for seed in ("neq_s1", "neq_s2"):
        w = W[seed]; tab[seed] = {}
        wa, wb = w[:h] / w[:h].sum(), w[h:] / w[h:].sum()
        for (a, b) in planes:
            A_, B_ = h2(Ap, a, b), h2(Bp, a, b); cnt = A_ + B_; T = cnt / cnt.sum(); mask = cnt < 2000
            E = h2(refp, a, b, w); E = E / E.sum()
            IA, IB = h2(refp[:h], a, b, wa), h2(refp[h:], a, b, wb)
            with np.errstate(all="ignore"):
                Z = np.where(mask, np.nan, np.log(E / T))
                N = np.where(mask, np.nan, np.log((A_ / A_.sum()) / (B_ / B_.sum())) / 2 + np.log((IA / IA.sum()) / (IB / IB.sum())) / 2)
            def wrms(Q): ok = ~np.isnan(Q); return float(np.sqrt(np.sum(T[ok] * Q[ok] ** 2) / np.sum(T[ok])))
            tab[seed][f"({lab[a]},{lab[b]})"] = (wrms(Z), wrms(N), int(mask.sum()), int(mask.size))
            if seed == "neq_s1" and (a, b) in rows:
                figrows.append((a, b, T, E, Z, N, cnt))
        print(n, seed, {k: (round(v[0], 3), round(v[1], 3)) for k, v in tab[seed].items()}, flush=True)
    R["marginals"] = tab
    heads = ["Simulation", "Neural network", "log(NN / simulation)", "noise reference"]
    fig, ax = plt.subplots(3, 4, figsize=(18, 12.5))
    for r, (a, b, T, E, Z, N, cnt) in enumerate(sorted(figrows, key=lambda t: rows.index((t[0], t[1])))):
        e = [*ext(a), *ext(b)]
        def wrms(Q): ok = ~np.isnan(Q); return np.sqrt(np.sum(T[ok] * Q[ok] ** 2) / np.sum(T[ok]))
        lt = np.log10(np.where(cnt > 0, T, np.nan)); vmin = np.nanpercentile(lt, 5); vmax = np.nanmax(lt)
        for k, Q in enumerate([lt, np.log10(np.where(E > 0, E, np.nan))]):
            im = ax[r, k].imshow(Q.T, origin="lower", extent=e, aspect="auto", cmap="jet", vmin=vmin, vmax=vmax)
            fig.colorbar(im, ax=ax[r, k], label="log₁₀ density"); ax[r, k].set_title(heads[k])
        for k, Q in enumerate([Z, N]):
            im = ax[r, k + 2].imshow(Q.T, origin="lower", extent=e, aspect="auto", cmap="coolwarm", vmin=-0.2, vmax=0.2)
            fig.colorbar(im, ax=ax[r, k + 2], label="log ratio"); ax[r, k + 2].set_title(f"{heads[k + 2]}   (rms {wrms(Q):.3f})")
        for k in range(4): ax[r, k].set_xlabel(lab[a]); ax[r, k].set_ylabel(lab[b])
    Tn = "T₃" if n == 3 else "T₄"
    fig.suptitle(f"2D marginals, (T₁,{Tn})=(2,8): neural density vs direct simulation", fontsize=14)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_marginals_n{n}_v5.png", dpi=110); plt.close(fig)
    # ---- angular slices at I_1 = ... = I_n = I, normalized on the plane (n = 4: theta_4 integrated out)
    m = M.load_bundle(f"{MOD}/n{n}_neq_s1")
    g = np.linspace(-np.pi, np.pi, 121); dA = (g[1] - g[0]) ** 2; sl = {}
    fig, ax = plt.subplots(2, 3, figsize=(19 if n == 3 else 16, 11.5 if n == 3 else 9.6))
    for c, I in enumerate((0.3, 0.5, 1.0)):
        if n == 3:
            G1, G3 = np.meshgrid(g, g, indexing="ij")
            xp = np.column_stack([np.full(G1.size, I)] * 3 + [G1.ravel(), G3.ravel()]); xf = xp.copy(); xf[:, 3] = -xf[:, 3]
            E = E_np(xf); fx = logratio(m, xf)
            dens = [(-E / 5 + fx).reshape(G1.shape), (-E / 5).reshape(G1.shape)]
            dens = [np.exp(L - L.max()) for L in dens]
        else:
            t4 = np.linspace(-np.pi, np.pi, 48, endpoint=False)
            G2, G3, G4 = np.meshgrid(g, g, t4, indexing="ij")
            xf = np.column_stack([np.full(G2.size, I)] * 4 + [G2.ravel(), G3.ravel(), G4.ravel()])
            E = E_np(xf); fx = logratio(m, xf)
            dens = []
            for L in ((-E / 5 + fx), (-E / 5)):
                L = L.reshape(G2.shape); dens.append(np.exp(L - L.max()).mean(-1))
        for r, p in enumerate(dens):
            p = p / (p.sum() * dA)
            im = ax[r, c].imshow(p.T, origin="lower", extent=[-np.pi, np.pi, -np.pi, np.pi], aspect="auto", cmap="jet")
            fig.colorbar(im, ax=ax[r, c]); ax[r, c].set_xlabel(lab[n]); ax[r, c].set_ylabel(lab[n + 1])
            ax[r, c].set_title((f"T₁=2, {Tn}=8 (NN)  I={I}" if r == 0 else f"T₁={Tn}=5 (Gibbs)  I={I}"))
            k = np.unravel_index(p.argmin(), p.shape); kk = np.unravel_index(p.argmax(), p.shape)
            sl[f"{'nn' if r == 0 else 'gibbs'}_I{I}"] = dict(min=(round(g[k[0]], 2), round(g[k[1]], 2)), max=(round(g[kk[0]], 2), round(g[kk[1]], 2)),
                                                           ratio=float(p.max() / p.min()), asym_exchange=float(np.abs(p - p.T).max() / p.max()),
                                                           asym_reflect=float(np.abs(p - p[::-1, ::-1]).max() / p.max()),
                                                           mass_both_neg=float(p[(g[:, None] < 0) & (g[None, :] < 0)].sum() * dA),
                                                           mass_both_pos=float(p[(g[:, None] > 0) & (g[None, :] > 0)].sum() * dA))
    fig.suptitle("Top: non-equilibrium (NN)  |  Bottom: equilibrium (Gibbs)", fontsize=14)
    fig.tight_layout(); fig.savefig(f"{OUT}/fig_theta_slices_n{n}_v5.png", dpi=110 if n == 3 else 120); plt.close(fig)
    R["slices"] = sl
    for k, v in sl.items(): print(n, k, v)
json.dump(res, open(f"{OUT}/an_sec4_v5.json", "w"), indent=1)
print("done")
