#!/usr/bin/env python3
"""Step 2 analysis: gEDMD spectra for nested dictionaries, half-sample cross-validation,
modal weights, and autocorrelation cross-check.  Usage: python3 analyze_step2.py PREFIX [PREFIX ...]

For every run prefix (OUT.bin/OUT.json from NLS3_gedmd_stream) writes
  PREFIX_analysis.json   eigenvalues per dictionary (full / half A / half B / per group),
                         cross-validated flags, weak-stationarity check, mean actions,
                         modal weights and predicted vs measured autocorrelations
  PREFIX_modes.npz       eigenvalues and coefficient vectors (top modes) per dictionary,
                         plus the D3-level matrices, for the eigenfunction plots of Step 3
"""
import json, sys
import numpy as np
import gedmd_n3 as G
import ged_io

DICTS = {"D1": (2, 1), "D2": (3, 2), "D3": (4, 2), "D4": (4, 3), "D5": (5, 2), "D6": (5, 3)}
NTOP = 12


def solve(A, B, cut=1e-10):
    return G.solve(A, B, cut)


def rq(c, A, B):
    """Rayleigh quotient of coefficient vector c on another sample's (A, B)."""
    return (np.conj(c) @ B @ c) / (np.conj(c) @ A @ c)


def obs_vectors(terms, idx, s):
    pos = {terms[k]: i for i, k in enumerate(idx)}
    K = len(idx)
    def v(pairs):
        a = np.zeros(K)
        for t, c in pairs:
            if t not in pos:
                return None
            a[pos[t]] += c
        return a
    z = (0, 0, 0)
    return {
        "I2": v([(((0, 1, 0), z, 0, 0, "c"), s[1])]),
        "I1+I3": v([(((1, 0, 0), z, 0, 0, "c"), s[0]), (((0, 0, 1), z, 0, 0, "c"), s[2])]),
        "I1-I3": v([(((1, 0, 0), z, 0, 0, "c"), s[0]), (((0, 0, 1), z, 0, 0, "c"), -s[2])]),
        "h1h2cos1": v([((z, (1, 1, 0), 1, 0, "c"), 1.0)]),
        "h2h3cos3": v([((z, (0, 1, 1), 0, 1, "c"), 1.0)]),
    }


def gibbs_mean(T, n=80, L=9.0):
    I = (np.arange(n) + 0.5) * L / n
    th = np.linspace(-np.pi, np.pi, 48, endpoint=False); c = np.cos(th)
    I1, I2, I3 = np.meshgrid(I, I, I, indexing="ij"); M = I1 + I2 + I3
    rho = np.exp(-(M ** 2 - 0.5 * (I1 ** 2 + I2 ** 2 + I3 ** 2)) / (2 * T))
    rho = rho * np.exp(-(2 * I1[..., None] * I2[..., None] * c) / (2 * T)).mean(-1)
    rho = rho * np.exp(-(2 * I2[..., None] * I3[..., None] * c) / (2 * T)).mean(-1)
    return [float((rho * v).sum() / rho.sum()) for v in (I1, I2, I3)]


def analyze(prefix):
    meta, terms, groups = ged_io.read(prefix)
    s = meta["scales"]; Gn = len(groups)
    full = ged_io.combine(groups, range(Gn))
    hA = ged_io.combine(groups, range(Gn // 2)); hB = ged_io.combine(groups, range(Gn // 2, Gn))
    per = [ged_io.combine(groups, [g]) for g in range(Gn)]
    res = {"prefix": prefix, "meta": {k: v for k, v in meta.items() if k != "terms"}, "n": full["n"],
           "mean_x": full["mean_x"].tolist(),
           "mean_I_by_group": [p["mean_x"][:3].tolist() for p in per]}
    if abs(meta["T1"] - meta["T3"]) < 1e-12:
        res["gibbs_mean_I"] = gibbs_mean(meta["T1"])
    # autocorrelation measured (all groups) and its spread across groups
    def acf_norm(c):
        if c["acfn"][0] == 0:  # --from runs carry no time ordering
            return np.full_like(c["acf"], np.nan), np.full(len(c["acf"]), np.nan)
        m = c["so"] / c["n"]; a = c["acf"] / c["acfn"][None] - (m ** 2)[:, None]
        return a / a[:, :1], a[:, 0]
    acf_full, var_full = acf_norm(full)
    acf_groups = np.array([acf_norm(p)[0] for p in per])
    lags = np.arange(meta["nlag"]) * meta["snapshot"]
    res["acf"] = {"lags": lags.tolist(), "observables": meta["observables"],
                  "measured": acf_full.tolist(), "group_se": (acf_groups.std(0) / np.sqrt(Gn)).tolist()}
    modes_out = {}
    res["dict"] = {}
    for name, (d, M) in DICTS.items():
        if d > meta["d"] or M > meta["M"]:
            continue
        idx = ged_io.sub_index(terms, d, M)
        sub = lambda c: (c["A"][np.ix_(idx, idx)], c["B"][np.ix_(idx, idx)])
        A, B = sub(full); AA, BA = sub(hA); AB, BB = sub(hB)
        lam, C, kept = solve(A, B)
        lamA, CA, _ = solve(AA, BA); lamB, CB, _ = solve(AB, BB)
        # cross-validation: modes from half A evaluated on half B and vice versa
        def xval(lam_h, C_h, Ao, Bo):
            out = []
            for j in range(min(NTOP, len(lam_h))):
                l2 = rq(C_h[:, j], Ao, Bo)
                out.append([float(lam_h[j].real), float(lam_h[j].imag), float(l2.real), float(l2.imag)])
            return out
        xa = xval(lamA, CA, AB, BB); xb = xval(lamB, CB, AA, BA)
        # full-sample modes are validated out of sample: the nearest mode of half A must reproduce
        # the eigenvalue when its eigenvector is evaluated on half B (and vice versa).
        def nearest(lam_h, z):
            return int(np.argmin(np.abs(lam_h - z)))
        valid = []
        for j in range(min(NTOP, len(lam))):
            tol = 0.05 * abs(lam[j]) + 0.05
            ja, jb = nearest(lamA, lam[j]), nearest(lamB, lam[j])
            qa = rq(CA[:, ja], AB, BB)   # half-A mode tested on half B
            qb = rq(CB[:, jb], AA, BA)   # half-B mode tested on half A
            ok = (lam[j].real <= 1e-6 and abs(lamA[ja] - lam[j]) < tol and abs(lamB[jb] - lam[j]) < tol
                  and abs(qa - lam[j]) < tol and abs(qb - lam[j]) < tol)
            valid.append(dict(re=float(lam[j].real), im=float(lam[j].imag),
                              halfA=[float(lamA[ja].real), float(lamA[ja].imag)], halfB=[float(lamB[jb].real), float(lamB[jb].imag)],
                              A_on_B=[float(qa.real), float(qa.imag)], B_on_A=[float(qb.real), float(qb.imag)],
                              validated=bool(ok)))
        per_group = []
        for p in per:
            lg, _, _ = solve(*sub(p))
            neg = [complex(v) for v in lg if v.real < -1e-6][:6]
            per_group.append([[v.real, v.imag] for v in neg])
        mL = full["mL"][idx]; sdL = full["sdL"][idx]
        entry = {"K": int(len(idx)), "kept": int(kept), "modes_full": valid,
                 "xval_A_on_B": xa, "xval_B_on_A": xb, "per_group_negative": per_group,
                 "weak_mean_over_sd_rms": float(np.sqrt(np.mean((mL[1:] / np.maximum(sdL[1:], 1e-300)) ** 2)))}
        # modal weights and predicted autocorrelation for observables inside the span
        ov = obs_vectors(terms, idx, s)
        entry["observables"] = {}
        good = [j for j, v in enumerate(valid) if v["validated"]]
        for oname, a in ov.items():
            if a is None:
                continue
            alpha = np.linalg.lstsq(C, a.astype(complex), rcond=None)[0]
            proj = (a @ A) @ C
            gbar = a @ A[:, 0]
            var = a @ A @ a - gbar ** 2
            pred = np.array([np.real(np.sum(alpha * np.exp(lam * t) * proj)) - gbar ** 2 for t in lags]) / var
            w = np.real(alpha * proj) / var
            oi = meta["observables"].index(oname)
            entry["observables"][oname] = {
                "weights_top": [float(w[j]) for j in range(min(NTOP, len(w)))],
                "predicted": pred.tolist(),
                "max_abs_diff": float(np.max(np.abs(pred - acf_full[oi])))}
        res["dict"][name] = entry
        modes_out[name + "_lam"] = lam[:NTOP]; modes_out[name + "_C"] = C[:, :NTOP]; modes_out[name + "_idx"] = idx
        if name == "D3":
            modes_out.update(D3_A=A, D3_B=B, D3_AA=AA, D3_BA=BA, D3_AB=AB, D3_BB=BB)
        first = [f"{v['re']:+.3f}{v['im']:+.2f}i{'' if v['validated'] else '(x)'}" for v in valid[1:6]]
        print(f"{prefix} {name} K={len(idx)}: {' '.join(first)}", flush=True)
    np.savez_compressed(prefix + "_modes.npz", scales=np.array(s), **modes_out)
    json.dump(res, open(prefix + "_analysis.json", "w"), indent=1)


if __name__ == "__main__":
    for p in sys.argv[1:]:
        analyze(p)
