# Analyse stationary indicator autocorrelations from local_test stat runs.
import json, numpy as np, sys
from scipy.optimize import curve_fit
def load(fn):
    d = json.load(open(fn)); NF = len(d["centers"]) * len(d["radii"]); NL = d["nlag"]
    S1 = np.array(d["S1"]); S2 = np.array(d["S2"]).reshape(len(S1), NF, NL); CN = np.array(d["CN"])
    npair = S1[:, NF]; hits = S1[:, :NF]
    return d, hits, S2, CN, npair
def corr(hits, S2, CN, npair):
    p = hits.sum(0) / npair.sum()                     # stationary probability
    Egg = S2.sum(0) / CN.sum(0)[None]                 # E[g(0) g(t)]
    C = (Egg - p[:, None] ** 2) / (p - p ** 2)[:, None]
    # batch-jackknife standard error
    B = len(hits); Cs = []
    for b in range(B):
        m = np.ones(B, bool); m[b] = False
        pb = hits[m].sum(0) / npair[m].sum(); Eb = S2[m].sum(0) / CN[m].sum(0)[None]
        Cs.append((Eb - pb[:, None] ** 2) / (pb - pb ** 2)[:, None])
    Cs = np.array(Cs); se = np.sqrt((B - 1) / B * ((Cs - Cs.mean(0)) ** 2).sum(0))
    return p, C, se
def model(t, a1, a2, c, d, s, w):
    return a1 * np.exp(-0.90 * t) + a2 * np.exp(-1.85 * t) + np.exp(-s * t) * (c * np.cos(w * t) + d * np.sin(w * t))
if __name__ == "__main__":
    res = {}
    for case in sys.argv[1:]:
        d, hits, S2, CN, npair = load(f"stat_{case}.json")
        p, C, se = corr(hits, S2, CN, npair); t = np.arange(d["nlag"]) * d["snapshot"]
        names = [f"c{c} r={r}" for c in range(len(d["centers"])) for r in d["radii"]]
        res[case] = dict(t=t.tolist(), names=names, centers=d["centers"], p=p.tolist(), C=C.tolist(), se=se.tolist(), fits={})
        print(f"== {case}: {npair.sum():.3g} (lane,snapshot) pairs")
        for f, nm in enumerate(names):
            c = C[f]; s = se[f]
            # turning points that are significant (difference > 3 se)
            tp = [round(t[i], 2) for i in range(2, len(t) - 2) if ((c[i] - c[i - 2]) * (c[i + 2] - c[i]) < 0) and abs(c[i] - c[i - 2]) > 3 * s[i] and abs(c[i + 2] - c[i]) > 3 * s[i]]
            m = (t >= 0.3) & (t <= 4.0)
            try:
                pf, pc = curve_fit(model, t[m], c[m], p0=[c[30], 0.0, 0.01, 0.0, 3.0, 3.5], sigma=s[m] + 1e-6, maxfev=20000)
                fit = dict(sigma=float(pf[4]), omega=float(abs(pf[5])), osc_amp=float(np.hypot(pf[2], pf[3])), slow_amp=float(pf[0]),
                           omega_se=float(np.sqrt(pc[5, 5])), sigma_se=float(np.sqrt(pc[4, 4])))
            except Exception as e:
                fit = dict(error=str(e))
            res[case]["fits"][nm] = fit
            print(f"{nm:10s} P={p[f]:.2e} hits={hits.sum(0)[f]:.0f}  C(0.1)={c[10]:.3f} C(0.5)={c[50]:.4f} C(1)={c[100]:.4f} C(2)={c[200]:.4f} "
                  f"min={c.min():+.4f}@{t[c.argmin()]:.2f}  turning pts={tp[:6]}  fit: {({k: round(v, 3) for k, v in fit.items()} if 'error' not in fit else fit)}")
    json.dump(res, open("stat_summary.json", "w"))
