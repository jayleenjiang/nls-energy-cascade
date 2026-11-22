"""Checks of frozen_phase_fp.solve: (1) I1 -> 0 limit vs the exact 1D stationary density; (2) grid refinement
N = 128, 256, 512 at three points.  (Monte Carlo cross-check: frozen_phase_sde.py, output in frozen_mc_check.txt.)"""
from frozen_phase_fp import solve, one_d_exact_sin
for model in ["canonical", "hlns"]:
    for I3 in [0.5, 2.0, 6.0]:
        ex = one_d_exact_sin(model, I3, 0.1, 8.0)
        r, _, _ = solve(model, 1e-4, I3, 0.1, 2.0, 8.0, N=256)
        print(f"{model} I3={I3}: exact 1D E[sin th3] = {ex:.6f}   2D FP with I1=1e-4: {r['E_sin3']:.6f}   diff {r['E_sin3'] - ex:+.1e}", flush=True)
for model in ["canonical", "hlns"]:
    for (I1, I3) in [(0.1, 6.0), (3.0, 3.0), (0.5, 2.0)]:
        vals = [solve(model, I1, I3, 0.1, 2.0, 8.0, N=N)[0]["F"] for N in (128, 256, 512)]
        print(f"{model} (I1,I3)=({I1},{I3}): F at N=128/256/512 = " + " / ".join(f"{v:.5f}" for v in vals), flush=True)
