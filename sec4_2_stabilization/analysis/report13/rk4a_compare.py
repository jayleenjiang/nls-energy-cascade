"""rk4a_chain (schemes4 RK4_A + Cartesian bath) vs aa_chain (Heun + action-angle EM bath) vs Cartesian v2, core 4.3 cases.
Usage: python3 rk4a_compare.py RK4A_RUNS AA_RUNS CART_RESULTS_JSON CART_DRAIN_JSON CART_TRAP_JSON OUTDIR"""
import json, sys
import numpy as np
RR, AA, CR, CD, CT, O = sys.argv[1:7]


def summ(p):
    d = {}
    for line in open(p):
        k, *v = line.split()
        try:
            d[k] = [float(x) for x in v]
        except ValueError:
            d[k] = v
    return d


def z(a, b):
    return (a[0] - b[0]) / np.hypot(a[1], b[1])


out = {"validation": {}, "case1": {}, "case3": {}}
s = summ(f"{RR}/val/n4_eqT5_full_summary.txt")
M = sum(s[f"meanI_{j}"][0] for j in range(1, 5)); Mse = np.sqrt(sum(s[f"meanI_{j}"][1] ** 2 for j in range(1, 5)))
out["validation"]["n4_full_T5_mean_M"] = {"rk4a": [M, Mse], "gibbs": 5.1246, "z": (M - 5.1246) / Mse}
s = summ(f"{RR}/val/n3_hold2_T5_5_summary.txt")
out["validation"]["n3_I2held_T5_meanI"] = {"I1": s["meanI_1"], "I3": s["meanI_3"], "exact": 1.7375,
                                           "z1": (s["meanI_1"][0] - 1.7375) / s["meanI_1"][1], "z3": (s["meanI_3"][0] - 1.7375) / s["meanI_3"][1]}
a, b = (summ(f"{RR}/val/n3_hold2_T2_8_{x}_summary.txt") for x in ("dt1e-3", "dt5e-4"))
out["validation"]["n3_T2_8_Lambda_dt"] = {"dt1e-3": a["rate_2"], "dt5e-4": b["rate_2"], "z": z(a["rate_2"], b["rate_2"])}
b4 = summ(f"{RR}/val/n4_hold2_T2_8_dt5e-4_summary.txt"); a4 = summ(f"{RR}/case1/n4_hold2_T2_8_summary.txt")
out["validation"]["n4_hold2_T2_8_Lambda_dt"] = {"dt1e-3": a4["rate_2"], "dt5e-4": b4["rate_2"], "z": z(a4["rate_2"], b4["rate_2"])}
nf = 0
cart = json.load(open(CR))["renorm"]
zc, za = [], []
for T in ("2_8", "5_5", "8_2"):
    for case, held, ck in (("n3_hold2", [2], "n3_z2"), ("n4_hold2", [2], "n4_z2"), ("n4_hold23", [2, 3], "n4_z23")):
        s = summ(f"{RR}/case1/{case}_T{T}_summary.txt"); sa = summ(f"{AA}/case1/{case}_T{T}_summary.txt")
        nf += s["nonfinite_steps"][0]
        cv = cart[f"{ck}_T{T}_dt1e-3"]
        for zz in held:
            lam, la, cl = s[f"rate_{zz}"], sa[f"rate_{zz}"], cv[f"site{zz}"]["Lambda_log"]
            r = {"Lambda_rk4a": lam, "Lambda_aa": la, "Lambda_cart": cl, "E_rk4a": [-lam[0] / 4, lam[1] / 4],
                 "z_rk4a_cart": z(lam, cl), "z_rk4a_aa": z(lam, la)}
            if case == "n4_hold23":
                r["E_single_rk4a"] = s["Sleft_2"] if zz == 2 else s["Sright_3"]
            zc.append(r["z_rk4a_cart"]); za.append(r["z_rk4a_aa"])
            out["case1"][f"{case}_T{T}_I{zz}"] = r
out["case1_summary"] = {"z_rms_vs_cart": float(np.sqrt(np.mean(np.square(zc)))), "max_abs_z_vs_cart": float(np.max(np.abs(zc))),
                        "z_rms_vs_aa": float(np.sqrt(np.mean(np.square(za)))), "max_abs_z_vs_aa": float(np.max(np.abs(za))),
                        "count": len(zc), "nonfinite_steps": nf}
s = summ(f"{RR}/ness/n3_T2_8_summary.txt")
Mn = sum(s[f"meanI_{j}"][0] for j in range(1, 4)); out["case3"]["ness_mean_M"] = Mn
cd = json.load(open(CD)); ct = json.load(open(CT))
for E0, dt in ((10, 1e-3), (30, 2e-4), (100, 1e-4), (200, 5e-5)):
    a = np.fromfile(f"{RR}/drain/n3_E{E0}_dump.f32", dtype=np.float32).reshape(16, -1, 16, 5).transpose(0, 2, 1, 3).reshape(256, -1, 5).astype(float)
    every = round(0.002 / dt); t = (np.arange(a.shape[1]) + 1) * every * dt
    I = a[:, :, :3]; Mt = I.sum(-1)
    th1, th3 = -a[:, :, 3], a[:, :, 4]
    tau = np.array([t[np.argmax(I[k, :, 1] <= E0 / 2)] if np.any(I[k, :, 1] <= E0 / 2) else np.nan for k in range(256)])
    mask = t[None, :] < tau[:, None]
    tauM = np.array([t[np.argmax(Mt[k] <= 2 * Mn)] if np.any(Mt[k] <= 2 * Mn) else np.inf for k in range(256)])
    r = {"tau_drain_median": float(np.nanmedian(tau)), "mean_sin_th1": float(np.mean(np.sin(th1[mask]))),
         "mean_sin_th3": float(np.mean(np.sin(th3[mask]))), "tau_return_median": float(np.median(tauM)),
         "not_returned_by_t10": int(np.isinf(tauM).sum())}
    k = f"n3_E{E0}"
    if k in cd:
        r["tau_drain_median_cart"] = cd[k]["tau_drain_median"]
    if k in ct:
        r["tau_return_median_cart"] = ct[k]["tau_M_median"]
    out["case3"][k] = r
print(json.dumps(out, indent=1, default=str))
json.dump(out, open(f"{O}/rk4a_compare.json", "w"), indent=1, default=str)
