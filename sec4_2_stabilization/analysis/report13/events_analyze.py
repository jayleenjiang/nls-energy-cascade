"""4.3.5: what happens in the actual NESS when a middle mode becomes small.
Inputs: events/{n3_T2_8,n3_T5_5,n4_T2_8}_dump.f32 (interval 0.01, 64 streams x 1000), ness/n{3,4}_T2_8_dump.f32
(interval 0.1, 256 streams x 2000).
(1) Small-action law: P(I_z < x) ~ x^kappa for small x (kappa = 1 exactly under the equal-temperature Gibbs law,
    where the density of I_z is finite and positive at 0).  kappa from a log-log fit of the empirical CDF on
    x in [1e-3, 3e-2]; SE from 16 groups of streams.
(2) Dips: a passage of I_z below b = 0.1 that reaches below a = 0.01.  For each dip: minimum, time from the
    minimum back to b, mean log-slope on the way up, and the bond angles (advisor-note convention) on the way down
    and on the way up.
Usage: python3 events_analyze.py HIDIR OUTDIR"""
import json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
H, O = sys.argv[1:3]


def load(path, n, nb):
    D = 2 * n - 1
    a = np.fromfile(path, dtype=np.float32).reshape(nb, -1, 16, D)
    return a.transpose(0, 2, 1, 3).reshape(nb * 16, -1, D).astype(np.float64)


def kappa(I, groups=16, lo=1e-3, hi=3e-2):
    xs = np.geomspace(lo, hi, 12)
    def fit(v):
        c = np.array([(v < x).mean() for x in xs])
        return np.polyfit(np.log(xs), np.log(c), 1)[0]
    k = fit(I.ravel())
    g = np.array([fit(part.ravel()) for part in np.array_split(I, groups, axis=0)])
    return float(k), float(g.std(ddof=1) / np.sqrt(groups)), [float((I < x).mean()) for x in xs], xs.tolist()


def note_angles(a, n, z):
    """angles on the two bonds of site z (one-based) in the advisor-note convention:
    left bond th_left = -theta_z (note th1 for z = 2), right bond th_right = theta_{z+1}"""
    left = -a[:, :, n + z - 2]
    right = a[:, :, n + z - 1] if z < n else None
    return left, right


def dips(I, left, right, Il, Ir, dt, a_lo=0.01, b_hi=0.1):
    ev = []; fall = {"l": [], "r": []}; rise = {"l": [], "r": []}
    for s in range(I.shape[0]):
        x = I[s]; below = x < b_hi
        k = 0; K = len(x)
        while k < K:
            if not below[k]:
                k += 1; continue
            k0 = k
            while k < K and below[k]:
                k += 1
            k1 = k                                          # [k0, k1) is a passage below b
            if k0 == 0 or k1 >= K:
                continue                                    # incomplete passage
            seg = x[k0:k1]
            if seg.min() >= a_lo:
                continue
            m = k0 + int(np.argmin(seg))
            up = (np.log(b_hi) - np.log(x[m])) / ((k1 - m) * dt) if k1 > m else np.nan
            ev.append({"min": float(x[m]), "t_down": (m - k0) * dt, "t_up": (k1 - m) * dt, "slope_up": float(up)})
            fall["l"].append(left[s, k0:m]); rise["l"].append(left[s, m:k1])
            if right is not None:
                fall["r"].append(right[s, k0:m]); rise["r"].append(right[s, m:k1])
    cat = lambda L: np.concatenate(L) if L else np.array([])
    return ev, {k: cat(v) for k, v in fall.items()}, {k: cat(v) for k, v in rise.items()}


res = {}
runs = [("n3_T2_8", 3, [2]), ("n3_T5_5", 3, [2]), ("n4_T2_8", 4, [2, 3])]
fig, axes = plt.subplots(2, 4, figsize=(17, 7.5))
col = 0
for name, n, zs in runs:
    a = load(f"{H}/events/{name}_dump.f32", n, 4)
    T_total = a.shape[0] * a.shape[1] * 0.01
    for z in zs:
        I = a[:, :, z - 1]
        k, kse, cdf, xs = kappa(I)
        left, right = note_angles(a, n, z)
        ev, fall, rise = dips(I, left, right, None, None, 0.01)
        mins = np.array([e["min"] for e in ev]); ups = np.array([e["slope_up"] for e in ev]); tup = np.array([e["t_up"] for e in ev])
        r = {"kappa_events": [k, kse], "dips": len(ev), "dips_per_unit_time_per_stream": len(ev) / T_total,
             "min_median": float(np.median(mins)) if len(ev) else None, "min_smallest": float(mins.min()) if len(ev) else None,
             "t_up_median": float(np.median(tup)) if len(ev) else None, "slope_up_median": float(np.median(ups)) if len(ev) else None,
             "mean_sin_left_fall": float(np.mean(np.sin(fall["l"]))), "mean_sin_left_rise": float(np.mean(np.sin(rise["l"])))}
        if n == 3 or z < n:
            r["mean_sin_right_fall"] = float(np.mean(np.sin(fall["r"]))); r["mean_sin_right_rise"] = float(np.mean(np.sin(rise["r"])))
        res[f"{name}_I{z}"] = r
        print(name, z, r, flush=True)
        if col < 4:
            ax = axes[0, col]
            for lab, arr, ls in [("left bond, falling", fall["l"], "--"), ("left bond, rising", rise["l"], "-"),
                                 ("right bond, falling", fall["r"], "--"), ("right bond, rising", rise["r"], "-")]:
                if len(arr) == 0:
                    continue
                h, e = np.histogram(arr, bins=72, range=(-np.pi, np.pi), density=True)
                ax.plot(0.5 * (e[1:] + e[:-1]), h, ls=ls, label=lab)
            for v in (2 * np.pi / 3, -2 * np.pi / 3):
                ax.axvline(v, color="grey", ls=":", lw=1)
            ax.set_title(f"{name}: angles while I{z} < 0.1 during a dip", fontsize=9); ax.set_xlabel("angle (notes convention)"); ax.legend(fontsize=6)
            ax = axes[1, col]
            ax.loglog(xs, cdf, "o-", ms=3, label=f"kappa = {k:.3f} ± {kse:.3f}")
            ax.loglog(xs, cdf[0] * (np.array(xs) / xs[0]), "k:", lw=1, label="slope 1 (Gibbs)")
            ax.set_xlabel(f"x"); ax.set_ylabel(f"P(I{z} < x)"); ax.set_title(f"{name}: small-I{z} law", fontsize=9); ax.legend(fontsize=7)
            col += 1
# long NESS runs (more independent samples) for kappa
for n, zs in [(3, [2]), (4, [2, 3])]:
    a = load(f"{H}/ness/n{n}_T2_8_dump.f32", n, 16)
    for z in zs:
        k, kse, _, _ = kappa(a[:, :, z - 1])
        res[f"ness_n{n}_T2_8_I{z}_kappa"] = [k, kse]
        print("ness", n, z, k, kse, flush=True)
    if n == 4:
        both = (a[:, :, 1] < 0.01) & (a[:, :, 2] < 0.01)
        res["ness_n4_T2_8_fraction_I2_and_I3_below_0.01"] = float(both.mean())
        res["ness_n4_T2_8_fraction_I2_below_0.01"] = float((a[:, :, 1] < 0.01).mean())
        res["ness_n4_T2_8_fraction_I3_below_0.01"] = float((a[:, :, 2] < 0.01).mean())
fig.tight_layout(); fig.savefig(f"{O}/fig_events.png", dpi=150)
json.dump(res, open(f"{O}/events.json", "w"), indent=1)
