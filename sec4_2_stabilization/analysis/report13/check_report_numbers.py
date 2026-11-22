"""Recompute every number printed in report_14 from the files of this package, section by section, and print
report value | value from the files | source.  Run from the package root:  python3 code/analysis/check_report_numbers.py
Reads only results/*.json, results/fokker_planck/*.json and data/ (summaries; dumps for the angle peaks, n=3 drain, NESS)."""
import json, numpy as np

def summ(p):
    d = {}
    for line in open(p):
        k, *v = line.split()
        try: d[k] = [float(x) for x in v]
        except ValueError: d[k] = v
    return d

rows = []
def row(sec, item, rep, val, src): rows.append((sec, item, rep, val, src))
L, U, R = "data/local_rk4a", "data/unity_batch06", "results"

# ---------------- n=3 case 1
S = "n=3 case 1"
for T in ("2_8", "5_5", "8_2"):
    f = f"{L}/case1/n3_hold2_T{T}_summary.txt"; lam, se = summ(f)["rate_2"]
    rep = {"2_8": ("-2.351", "9.40"), "5_5": ("-2.153", "8.61"), "8_2": ("-2.357", "9.43")}[T]
    row(S, f"E, ({T.replace('_', ',')})", rep[0], f"{-lam / 4:.3f} ± {se / 4:.3f}", f"{f} → rate_2 (E = −rate_2/4)")
    row(S, f"Λ, ({T.replace('_', ',')})", rep[1], f"{lam:.3f} ± {se:.3f}", f"{f} → rate_2")
Ts = ["0.5", "1", "2", "5", "10"]
lam = [summ(f"{L}/sweep/n3_hold2_T{a}_{b}_summary.txt")["rate_2"][0] for a in Ts for b in Ts]
row(S, "sweep, 25 temperature pairs", "Λ = 2.7–12.1, all > 0", f"{min(lam):.3f}–{max(lam):.3f}", f"{L}/sweep/*_summary.txt → rate_2 (also {R}/sweep_rk4a.json)")
modes = []
for T in ("2_8", "5_5", "8_2"):
    a = np.fromfile(f"{L}/angles/n3_hold2_T{T}_dump.f32", dtype=np.float32).reshape(-1, 5).astype(float)
    for th in (-a[:, 3], a[:, 4]):
        h, e = np.histogram(np.remainder(th + np.pi, 2 * np.pi) - np.pi, bins=360, range=(-np.pi, np.pi))
        modes.append(0.5 * (e[1:] + e[:-1])[np.convolve(h, np.ones(9) / 9, "same").argmax()])
row(S, "angle peak (θ1, θ3, three temperatures)", "4π/3 ≡ −2.094", f"{min(modes):.3f} to {max(modes):.3f}", f"{L}/angles/n3_hold2_*_dump.f32 (figure: {R}/fig43_boundary_angles_rk4a.png)")

# ---------------- n=3 case 2
S = "n=3 case 2"
rep = {"canonical_T2_8": ["-0.15", "-0.43", "-0.70", "-0.85"], "canonical_T5_5": ["-0.23", "-0.53", "-0.77", "-0.86"], "hlns": ["-0.83", "-0.85", "-0.86", "-0.86"]}
for k in rep:
    d = json.load(open(f"{R}/fokker_planck/frozen_map_part_{k}.json")); g = np.array(d["grid"]); F = np.array(d["maps"][k]["F"])
    for r_, I3 in zip(rep[k], (0.5, 1.0, 2.0, 6.0)):
        j = int(np.argmin(abs(g - I3)))
        row(S, f"F/I3, {k}, I1=0.05, I3={I3:g}", r_, f"{F[0, j] / I3:.3f}", f"{R}/fokker_planck/frozen_map_part_{k}.json → maps.{k}.F[0][{j}] / I3")

# ---------------- n=3 case 3
S = "n=3 case 3"
s = summ(f"{L}/ness/n3_T2_8_summary.txt"); Mn = sum(s[f"meanI_{j}"][0] for j in (1, 2, 3))
row(S, "<M>_NESS (n=3)", "4.70", f"{Mn:.3f}", f"{L}/ness/n3_T2_8_summary.txt → meanI_1+meanI_2+meanI_3")
rep = {10: ("0.168", "0.50", "0.44", "0.21"), 30: ("0.076", "0.59", "0.57", "1.11"), 100: ("0.030", "0.70", "0.68", "1.65"), 200: ("0.016", "0.74", "0.72", "1.78")}
for E0, dt in ((10, 1e-3), (30, 2e-4), (100, 1e-4), (200, 5e-5)):
    f = f"{L}/drain/n3_E{E0}_dump.f32"
    a = np.fromfile(f, dtype=np.float32).reshape(16, -1, 16, 5).transpose(0, 2, 1, 3).reshape(256, -1, 5).astype(float)
    t = (np.arange(a.shape[1]) + 1) * round(0.002 / dt) * dt; I = a[:, :, :3]; M = I.sum(-1)
    tau = np.array([t[np.argmax(x <= E0 / 2)] for x in I[:, :, 1]]); mask = t[None, :] < tau[:, None]
    tauM = np.array([t[np.argmax(m <= 2 * Mn)] if (m <= 2 * Mn).any() else np.inf for m in M])
    v = (f"{np.median(tau):.3f}", f"{np.mean(np.sin(-a[:, :, 3][mask])):.2f}", f"{np.mean(np.sin(a[:, :, 4][mask])):.2f}", f"{np.median(tauM):.2f}")
    for nm, rr, vv in zip(("first t with I2<=E0/2", "mean sin θ1", "mean sin θ3", "first t with M<=2<M>"), rep[E0], v):
        row(S, f"E0={E0}: {nm}", rr, vv, f"{f} (256 trajectories; times are medians, sin are means)")

# ---------------- n=4 cases 1, 2
S = "n=4 case 1"
rep = {"2_8": ("-1.525", "6.10"), "5_5": ("-1.568", "6.27"), "8_2": ("-2.024", "8.10")}
for T in rep:
    f = f"{L}/case1/n4_hold2_T{T}_summary.txt"; lam, se = summ(f)["rate_2"]
    row(S, f"E, ({T.replace('_', ',')})", rep[T][0], f"{-lam / 4:.3f}", f"{f} → rate_2"); row(S, f"Λ2, ({T.replace('_', ',')})", rep[T][1], f"{lam:.3f} ± {se:.3f}", f"{f} → rate_2")
S = "n=4 case 2"
rep = {"2_8": ("-0.422", "-2.135", "1.69 / 8.54"), "5_5": ("-1.241", "-1.243", "4.97 / 4.97"), "8_2": ("-2.138", "-0.423", "8.55 / 1.69")}
for T in rep:
    f = f"{L}/case1/n4_hold23_T{T}_summary.txt"; s = summ(f)
    row(S, f"E{{I1 sin θ1}}, ({T.replace('_', ',')})", rep[T][0], f"{s['Sleft_2'][0]:.3f}", f"{f} → Sleft_2")
    row(S, f"E{{I4 sin θ4}}, ({T.replace('_', ',')})", rep[T][1], f"{s['Sright_3'][0]:.3f}", f"{f} → Sright_3")
    row(S, f"Λ2 / Λ3, ({T.replace('_', ',')})", rep[T][2], f"{s['rate_2'][0]:.3f} / {s['rate_3'][0]:.3f}", f"{f} → rate_2, rate_3")

# ---------------- n=4 case 3
S = "n=4 case 3"; ts = json.load(open(f"{R}/trap_scan.json"))
rep = {50: ("0.48", "11", "-0.37"), 70: ("0.67", "21", "-0.63"), 100: ("0.72", "39", "-0.76"), 140: ("0.66", "68", "-0.81"), 200: ("0.66", "116", "-0.78"), 280: ("0.69", "183", "-0.79"), 400: ("0.67", "272", "-0.78")}
for E0, (r1, r2, r3) in rep.items():
    e = ts[f"e0scan/n4_E{E0}"]; src = f"{R}/trap_scan.json → e0scan/n4_E{E0}"
    row(S, f"E0={E0}: fraction at t=5", r1, f"{e['formed_fraction_t5']:.2f}", src + ".formed_fraction_t5")
    row(S, f"E0={E0}: median first t with M<=2<M>", r2, f"{e['tau_M_median_formed']:.1f}", src + ".tau_M_median_formed")
    row(S, f"E0={E0}: mean dM/dt", r3, f"{e['dMdt_fd_pair']:.3f} ± {e['dMdt_fd_pair_se']:.3f}", src + ".dMdt_fd_pair")
big = [ts[f"e0scan/n4_E{E}"] for E in (70, 100, 140, 200, 280, 400)]
row(S, "fraction range, E0>=70", "66–72%", f"{min(b['formed_fraction_t5'] for b in big):.2f}–{max(b['formed_fraction_t5'] for b in big):.2f}", f"{R}/trap_scan.json → formed_fraction_t5")
row(S, "M at t=5 / E0", "≈0.57", f"{min(b['M_at_t5_formed_mean'] / b['E0'] for b in big):.3f}–{max(b['M_at_t5_formed_mean'] / b['E0'] for b in big):.3f}", f"{R}/trap_scan.json → M_at_t5_formed_mean / E0")
row(S, "M when the state ends", "≈16–17", f"{min(b['M_at_pair_end_median'] for b in big):.1f}–{max(b['M_at_pair_end_median'] for b in big):.1f}", f"{R}/trap_scan.json → M_at_pair_end_median")
row(S, "<M>_NESS (n=4)", "5.24", f"{big[0]['M_ref']:.3f}", f"{R}/trap_scan.json → M_ref (= mean M of {U}/ness/n4_T2_8)")
tl = json.load(open(f"{R}/tail.json"))["n4_NESS"]
row(S, "NESS: largest M", "24.7 in 5.12e6 samples", f"{tl['max']:.2f} in {tl['samples']}", f"{R}/tail.json → n4_NESS.max, samples")
tr = json.load(open(f"{R}/trap.json"))
row(S, "long stay needs E0 ≳ 50", "E0 ≳ 50", ", ".join(f"E0={E}: {tr[f'n4_E{E}']['tau_M_median']:.1f}" for E in (20, 30, 50, 100, 200)), f"{R}/trap.json → n4_E*.tau_M_median (64 traj., {U}/trap)")

# ---------------- n=5
S = "n=5"; n5 = json.load(open(f"{R}/n5.json"))
runs = [k for k in n5 if k.startswith("n5_I")]
l23 = [n5[k]["leak_by_pair"]["2-3"]["dMdt_fd"] for k in runs]
l34 = [n5[k]["leak_by_pair"]["3-4"]["dMdt_fd"] for k in runs if not k.endswith("E50")]
row(S, "<M>_NESS (n=5)", "5.73", f"{n5['n5_ness_mean_M']:.3f}", f"{R}/n5.json → n5_ness_mean_M")
row(S, "mean dM/dt, I2,I3 (six runs)", "−0.161 to −0.189", f"{max(l23):.3f} to {min(l23):.3f}", f"{R}/n5.json → n5_I*_E*.leak_by_pair.2-3.dMdt_fd")
row(S, "mean dM/dt, I3,I4 (E0>=100)", "−0.660 to −0.713", f"{max(l34):.3f} to {min(l34):.3f}", f"{R}/n5.json → n5_I*_E100/200.leak_by_pair.3-4.dMdt_fd")
t23 = [n5[f"n5_{w}_E100"]["cluster_at_t5"]["2-3"]["tau_M_median"] for w in ("I2", "I3")]
t34 = [n5[f"n5_{w}_E100"]["cluster_at_t5"]["3-4"]["tau_M_median"] for w in ("I2", "I3")]
row(S, "median first t with M<=2<M>, I2,I3, E0=100", "269–344", f"{min(t23):.1f}–{max(t23):.1f}", f"{R}/n5.json → n5_I*_E100.cluster_at_t5.2-3.tau_M_median")
row(S, "median first t with M<=2<M>, I3,I4, E0=100", "70–71", f"{min(t34):.1f}–{max(t34):.1f}", f"{R}/n5.json → n5_I*_E100.cluster_at_t5.3-4.tau_M_median")
pairs = [n5[k]["cluster_while_M>20_after_t5"].get("2-3", {}).get("time_fraction", 0) + n5[k]["cluster_while_M>20_after_t5"].get("3-4", {}).get("time_fraction", 0) for k in runs]
row(S, "time with M>20 in I2,I3 or I3,I4", "≥ 98.6%", f"min {min(pairs):.3f}", f"{R}/n5.json → cluster_while_M>20_after_t5.(2-3, 3-4).time_fraction")

sec = None
for s_, it, rp, v, src in rows:
    if s_ != sec:
        print(f"\n### {s_}\n\n| item | report | from files | source |\n|---|---|---|---|"); sec = s_
    print(f"| {it} | {rp} | {v} | `{src}` |")
