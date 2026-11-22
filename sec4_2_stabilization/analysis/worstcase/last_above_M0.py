"""Last t with M(t) > M(0): (a) for the mean curve E M(t); (b) per trajectory (median, 90% quantile, max).
Data: run_Mrise.sh (n = 3..6, gamma = 1e-3 / 1e-4, t <= 0.1, sampled every 1e-3). Writes last_above_M0.json."""
import json
import numpy as np
iv = 1e-3; M0 = 100.02; res = {}
for g in ("1e-3", "1e-4"):
    for n in (3, 4, 5, 6):
        for ph in ("worst", "random"):
            a = np.fromfile(f"n{n}_g{g}_{ph}_dump.f32", dtype=np.float32).reshape(-1, 100, 16, 2 * n - 1)
            M = a[..., :n].astype(float).sum(-1).transpose(1, 0, 2).reshape(100, -1)      # (time, trajectories)
            t = (np.arange(100) + 1) * iv
            m = M.mean(1); above = np.where(m > M0)[0]
            mean_last = float(t[above[-1]]) if len(above) else 0.0
            ab = M > M0; anyab = ab.any(0); last = np.where(anyab, t[99 - np.argmax(ab[::-1], 0)], 0.0)
            res[f"n{n}_g{g}_{ph}"] = {"trajectories": int(M.shape[1]), "last_t_mean_curve": mean_last,
                                      "per_traj_median": float(np.median(last)), "per_traj_q90": float(np.quantile(last, .9, method="nearest")),
                                      "per_traj_max": float(last.max()), "frac_never_above": float((~anyab).mean()),
                                      "frac_above_at_end": float(ab[-1].mean())}
            print(n, g, ph, res[f"n{n}_g{g}_{ph}"], flush=True)
json.dump(res, open("last_above_M0.json", "w"), indent=1)
