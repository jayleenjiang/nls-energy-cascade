"""Numbers of report_13, n = 3 cases 1 and 3 and n = 4 cases 1 and 2, from the local rk4a_chain runs only.
Usage: python3 tables_n3_n4.py LOCAL_RK4A_DIR   (folder with case1/, sweep/, drain/, ness/)
  case 1:  Lambda = time average of d log I_z/dt of the held mode (summary key rate_z), E = -Lambda/4, SE over streams.
  n = 4, I2 = I3 = 0:  E{I1 sin th1} = Sleft_2, E{I4 sin th4} = Sright_3 (notes convention).
  case 3 (n = 3 drain): drain time = first t with I2 <= E0/2 (median); mean sin th1, sin th3 while I2 > E0/2;
  return time = first t with M <= 2 <M>_NESS (median), <M>_NESS from ness/n3_T2_8."""
import sys, numpy as np
R = sys.argv[1]


def summ(p):
    d = {}
    for line in open(p):
        k, *v = line.split()
        try:
            d[k] = [float(x) for x in v]
        except ValueError:
            d[k] = v
    return d


for case, held in (("n3_hold2", [2]), ("n4_hold2", [2]), ("n4_hold23", [2, 3])):
    for T in ("2_8", "5_5", "8_2"):
        s = summ(f"{R}/case1/{case}_T{T}_summary.txt")
        for z in held:
            L, se = s[f"rate_{z}"]
            line = f"{case} ({T.replace('_', ',')}) I{z}: Lambda = {L:.3f} +- {se:.3f}   E = -Lambda/4 = {-L / 4:.3f} +- {se / 4:.3f}"
            if case == "n4_hold23":
                v = s["Sleft_2"] if z == 2 else s["Sright_3"]
                line += f"   single term ({'E[I1 sin th1]' if z == 2 else 'E[I4 sin th4]'}) = {v[0]:.3f} +- {v[1]:.3f}"
            print(line)
Ts = ["0.5", "1", "2", "5", "10"]
L = [summ(f"{R}/sweep/n3_hold2_T{a}_{b}_summary.txt")["rate_2"][0] for a in Ts for b in Ts]
print(f"n3 temperature sweep (25 runs): Lambda from {min(L):.3f} to {max(L):.3f}, all positive: {min(L) > 0}")
s = summ(f"{R}/ness/n3_T2_8_summary.txt")
Mn = sum(s[f"meanI_{j}"][0] for j in (1, 2, 3)); Mse = np.sqrt(sum(s[f"meanI_{j}"][1] ** 2 for j in (1, 2, 3)))
print(f"n3 <M>_NESS (2,8) = {Mn:.3f} +- {Mse:.3f}")
for E0, dt in ((10, 1e-3), (30, 2e-4), (100, 1e-4), (200, 5e-5)):
    a = np.fromfile(f"{R}/drain/n3_E{E0}_dump.f32", dtype=np.float32).reshape(16, -1, 16, 5).transpose(0, 2, 1, 3).reshape(256, -1, 5).astype(float)
    every = round(0.002 / dt); t = (np.arange(a.shape[1]) + 1) * every * dt
    I = a[:, :, :3]; M = I.sum(-1); th1, th3 = -a[:, :, 3], a[:, :, 4]
    tau = np.array([t[np.argmax(I[k, :, 1] <= E0 / 2)] if np.any(I[k, :, 1] <= E0 / 2) else np.nan for k in range(256)])
    mask = t[None, :] < tau[:, None]
    tauM = np.array([t[np.argmax(M[k] <= 2 * Mn)] if np.any(M[k] <= 2 * Mn) else np.inf for k in range(256)])
    print(f"n3 drain E0={E0}: drain time {np.nanmedian(tau):.3f}, mean sin th1 {np.mean(np.sin(th1[mask])):.2f}, "
          f"mean sin th3 {np.mean(np.sin(th3[mask])):.2f}, return time {np.median(tauM):.2f}, not returned by t=10: {int(np.isinf(tauM).sum())}")
