import json, numpy as np, sys
sys.path.insert(0, '..'); sys.path.insert(0, '.')
import gedmd_n3 as G
R = None   # folder holding CASE_modes.npz and CASE.json (with trailing /); set by the caller, e.g. E.R = "step2_run/"

def load_modes(case, D="D6"):
    Z = np.load(R + case + "_modes.npz", allow_pickle=True)
    meta = json.load(open(R + case.split("_dt")[0] + "_dt1e-3.json"))
    terms_full = [((t[0], t[1], t[2]), (t[3], t[4], t[5]), t[6], t[7], "s" if t[8] else "c") for t in meta["terms"]]
    idx = Z[D + "_idx"]; terms = [terms_full[i] for i in idx]
    return Z[D + "_lam"], Z[D + "_C"], terms, Z["scales"]

def phi_eval(x, terms, s, chunk=20000):
    out = []
    for i in range(0, len(x), chunk):
        P, *_ = G.evaluate(x[i:i + chunk], terms, s, 2.0, 8.0)   # Phi does not depend on T
        out.append(P)
    return np.vstack(out)

def real_modes(lam, C, Phi_samples, I2):
    """psi_j on samples; fix complex phase, normalise to unit stationary variance, sign: corr(psi, I2) > 0."""
    psis = []
    for j in range(C.shape[1]):
        v = Phi_samples @ C[:, j]
        ph = np.exp(-1j * np.angle(np.sum(v * v)) / 2)       # make it as real as possible
        c = C[:, j] * ph
        v = (Phi_samples @ c)
        sc = np.sqrt(np.mean(np.abs(v - v.mean()) ** 2))
        c = c / sc
        v = Phi_samples @ c
        if np.real(np.mean((v - v.mean()) * (I2 - I2.mean()))) < 0: c = -c
        psis.append(c)
    return np.array(psis).T

def psi_eval(x, terms, s, C, chunk=4000):
    out = []
    for i in range(0, len(x), chunk):
        P, *_ = G.evaluate(x[i:i + chunk], terms, s, 2.0, 8.0)
        out.append(P @ C)
    return np.vstack(out)

def realify(C, V, I2):
    """Fix phase, unit stationary variance, sign corr(psi, I2) >= 0. Returns (C_r, V_r)."""
    C = C.copy().astype(complex); V = V.copy().astype(complex)
    for j in range(C.shape[1]):
        v = V[:, j] - V[:, j].mean()
        ph = np.exp(-1j * np.angle(np.sum(v * v)) / 2)
        sc = np.sqrt(np.mean(np.abs(v) ** 2))
        f = ph / sc
        if np.real(np.mean(v * f * (I2 - I2.mean()))) < 0: f = -f
        C[:, j] *= f; V[:, j] *= f
    return C, V
