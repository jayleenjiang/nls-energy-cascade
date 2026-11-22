"""4.3 extension, B2: how often a high-M state of the n = 4 chain is an anti-phase pair, NESS (2,8) vs Gibbs.
At fixed M, H = M^2 - 1/2 sum I^2 + 2 sum I_{j-1} I_j cos th_j is minimised by an anti-phase pair on any bond
(H = M^2 - 3M^2/4), so Gibbs exp(-H/2T) at large M concentrates on pairs; the bath only acts on the endpoints.
Gibbs: Metropolis in (log I_1..4, th_2..4) (phases uniform => bond angles independent uniform a priori).
Only states with M > 10 are stored. Pair on bond b (modes b-1, b): (I_{b-1} + I_b)/M > 0.9 and cos th_b < -0.9.
Usage: python3 pair_gibbs.py HIDIR OUTDIR"""
import json, sys
import numpy as np
H_, O = sys.argv[1:3]


def ham(I, th):
    M = I.sum(-1)
    return M ** 2 - 0.5 * (I ** 2).sum(-1) + 2 * (I[..., :-1] * I[..., 1:] * np.cos(th)).sum(-1)


def classify(I, th, Mmin):
    M = I.sum(-1); sel = M > Mmin
    out = {"n_states": int(sel.sum())}
    for b, name in [(0, "pair_12"), (1, "pair_23"), (2, "pair_34")]:
        p = ((I[:, b] + I[:, b + 1]) / M > 0.9) & (np.cos(th[:, b]) < -0.9)
        out[name] = float(p[sel].mean()) if sel.any() else None
    return out


def gibbs(T, chains=4096, burn=3000, steps=30000, thin=5, seed=7):
    rng = np.random.default_rng(seed)
    x = np.log(np.full((chains, 4), np.sqrt(T) * 0.5)); th = rng.uniform(-np.pi, np.pi, (chains, 3))
    lp = -ham(np.exp(x), th) / (2 * T) + x.sum(-1)
    I_out, th_out = [], []; msum = 0.0; mcount = 0
    for k in range(burn + steps):
        y = x + 0.3 * rng.standard_normal(x.shape); ph = th + 0.6 * rng.standard_normal(th.shape)
        lq = -ham(np.exp(y), ph) / (2 * T) + y.sum(-1)
        acc = np.log(rng.random(chains)) < lq - lp
        x[acc] = y[acc]; th[acc] = ph[acc]; lp[acc] = lq[acc]
        if k >= burn and k % thin == 0:
            E = np.exp(x); M = E.sum(-1); msum += M.sum(); mcount += chains
            keep = M > 10                     # store only the high-M states (memory)
            I_out.append(E[keep]); th_out.append(th[keep].copy())
    I = np.concatenate(I_out); th = (np.concatenate(th_out) + np.pi) % (2 * np.pi) - np.pi
    return I, th, msum / mcount, mcount


res = {}
a = np.fromfile(f"{H_}/ness/n4_T2_8_dump.f32", dtype=np.float32).reshape(-1, 7).astype(float)
I, th = a[:, :4], a[:, 4:]                          # th_2, th_3, th_4 (long-chain convention)
for Mmin in (10, 15, 20):
    res[f"NESS_2_8_M>{Mmin}"] = classify(I, th, Mmin)
for T in (5.0, 8.0):
    I, th, mM, ns = gibbs(T)
    res[f"Gibbs_T{T:g}_mean_M"] = float(mM); res[f"Gibbs_T{T:g}_samples"] = int(ns)
    for Mmin in (10, 15, 20):
        res[f"Gibbs_T{T:g}_M>{Mmin}"] = classify(I, th, Mmin)
    print(T, {k: v for k, v in res.items() if f"T{T:g}" in k}, flush=True)
print(json.dumps(res, indent=1))
json.dump(res, open(f"{O}/pair_gibbs.json", "w"), indent=1)
