"""4.3.3 maps: F(I1, I3) for canonical bath ((T1,T3) = (2,8) and (5,5), gamma = 0.1) and HLNS reduced dynamics
(gamma = 0.1), N = 512 (with N = 256 recorded as a resolution estimate), plus the averaging check: mean of
F_canonical over (I1, I3) drawn from the restricted Gibbs law at T = 5 (400 draws, N = 256).
One part per call so parts can run in parallel (and be resumed):
Usage: python3 frozen_phase_map.py OUT.json {canonical_T2_8 | canonical_T5_5 | hlns | averaging}"""
import json, sys
import numpy as np
from frozen_phase_fp import solve
grid = [0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 3.0, 4.5, 6.0]
cases = {"canonical_T2_8": ("canonical", 2.0, 8.0), "canonical_T5_5": ("canonical", 5.0, 5.0), "hlns": ("hlns", 2.0, 8.0)}
part = sys.argv[2]
out = {"grid": grid, "maps": {}}
if part in cases:
    model, T1, T3 = cases[part]
    F = np.zeros((len(grid), len(grid))); F256 = F.copy(); S1 = F.copy(); S3 = F.copy()
    for i, I1 in enumerate(grid):
        for j, I3 in enumerate(grid):
            r, _, _ = solve(model, I1, I3, 0.1, T1, T3, N=512)
            r2, _, _ = solve(model, I1, I3, 0.1, T1, T3, N=256)
            F[i, j], F256[i, j], S1[i, j], S3[i, j] = r["F"], r2["F"], r["E_sin1"], r["E_sin3"]
        print(part, I1, np.round(F[i], 3).tolist(), flush=True)
    out["maps"][part] = {"F": F.tolist(), "F_N256": F256.tolist(), "E_sin1": S1.tolist(), "E_sin3": S3.tolist()}
elif part == "averaging":
    # draw (I1, I3) from exp(-[(I1+I3)^2 - (I1^2+I3^2)/2] / (2T)), T = 5, on [0, inf)^2
    rng = np.random.default_rng(7)
    h = 0.005; I = np.arange(h / 2, 15, h)
    A, B = np.meshgrid(I, I, indexing="ij")
    w = np.exp(-((A + B) ** 2 - 0.5 * (A ** 2 + B ** 2)) / 10.0).ravel(); w /= w.sum()
    pick = rng.choice(w.size, size=400, p=w)
    vals = []
    for c, k in enumerate(pick):
        I1, I3 = A.ravel()[k] + rng.uniform(-h / 2, h / 2), B.ravel()[k] + rng.uniform(-h / 2, h / 2)
        r, _, _ = solve("canonical", I1, I3, 0.1, 5.0, 5.0, N=256)
        vals.append(r["F"])
        if c % 50 == 0:
            print("averaging", c, flush=True)
    vals = np.array(vals)
    out["averaging_T5"] = {"mean_F": float(vals.mean()), "se": float(vals.std(ddof=1) / np.sqrt(len(vals))), "draws": len(vals)}
    print("averaging check T=5:", out["averaging_T5"], flush=True)
json.dump(out, open(sys.argv[1], "w"), indent=1)
