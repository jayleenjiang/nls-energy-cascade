"""Visualisation of n = 4, case 3 (I2 = E0 start): one trapped trajectory, the action profile, and M(t).
Data: Unity batch_06 e0scan (rk4a_chain), n = 4 NESS run, and the local n = 3 drain run (E0 = 200)."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
U, L, OUT = sys.argv[1:4]           # Unity out/, local rk4a runs/, output png
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]   # categorical slots 1-4 (validated)
INK, MUTED, GRID = "#2b2b2b", "#6b6b6b", "#e6e6e3"
plt.rcParams.update({"font.size": 9, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "xtick.color": MUTED,
                     "ytick.color": MUTED, "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.6})

def load(path, n, nstreams, iv):
    D = 2 * n - 1
    a = np.fromfile(path, dtype=np.float32).reshape(nstreams // 16, -1, 16, D).transpose(0, 2, 1, 3).reshape(nstreams, -1, D).astype(float)
    return a, (np.arange(a.shape[1]) + 1) * iv

ness = np.fromfile(f"{U}/ness/n4_T2_8_dump.f32", dtype=np.float32).reshape(-1, 7).astype(float)
Mref = ness[:, :4].sum(-1).mean(); prof_ness = ness[:, :4].mean(0)

def trap_stats(E0):
    a, t = load(f"{U}/e0scan/n4_E{E0}_dump.f32", 4, 128, 0.05)
    I = a[:, :, :4]; M = I.sum(-1); th3 = a[:, :, 5]
    pair = ((I[:, :, 1] + I[:, :, 2]) / M > 0.9) & (np.cos(th3) < -0.9)
    formed = pair[:, 99]
    tau = np.array([t[np.argmax(m <= 2 * Mref)] if (m <= 2 * Mref).any() else np.inf for m in M])
    return a, t, I, M, th3, pair, formed, tau

fig = plt.figure(figsize=(12.5, 8.2))
gs = fig.add_gridspec(3, 2, height_ratios=[2.3, 0.9, 0.15], width_ratios=[1.25, 1], hspace=0.45, wspace=0.22)
axA = fig.add_subplot(gs[0, 0]); axB = fig.add_subplot(gs[1, 0], sharex=axA)
gsr = gs[0:2, 1].subgridspec(2, 1, hspace=0.45); axC = fig.add_subplot(gsr[0]); axD = fig.add_subplot(gsr[1])

# ---- A, B: one trapped trajectory, E0 = 200 (lifetime closest to the median of the trapped ones)
a, t, I, M, th3, pair, formed, tau = trap_stats(200)
c = np.where(formed)[0]; s = c[np.argmin(np.abs(tau[c] - np.median(tau[c])))]
T0 = np.r_[0.0, t]; Is = np.vstack([[0.01, 200, 0.01, 0.01], I[s]]); Ms = Is.sum(-1)
xmax = 140
for j in range(4):
    axA.plot(T0, Is[:, j], color=C[j], lw=1.3 if j in (1, 2) else 1.6, zorder=3)
axA.plot(T0, Ms, color=INK, ls="--", lw=1.4, zorder=4)
axA.set_xlim(0, xmax); axA.set_ylim(-3, 205)
axA.text(72, Ms[int(70 / 0.05)] + 9, "M (total)", color=INK, fontsize=9, fontweight="bold")
axA.text(20, 75, r"$I_2$", color=C[1], fontsize=10, fontweight="bold"); axA.text(26, 75, r"$I_3$", color=C[2], fontsize=10, fontweight="bold")
axA.text(20, 7, r"$I_1$, $I_4$ (endpoints) stay $\approx 0$", color=MUTED, fontsize=8.5, va="bottom")
tb = tau[s]
axA.axvline(tb, color=MUTED, lw=0.8, ls=":")
axA.annotate("start: all action in $I_2$ ($E_0=200$).\nWithin $t<0.5$: ~37% leaves through the\nendpoints, the rest splits equally into $I_2$, $I_3$",
             xy=(0.3, 150), xytext=(10, 178), fontsize=8.5, color=INK, arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
axA.text(35, 140, "anti-phase pair $I_2\\approx I_3$:\nM decreases slowly and steadily", fontsize=8.5, color=INK)
axA.text(tb + 1.5, 140, f"pair breaks up\n(t ≈ {tb:.0f}); the rest\ndrains quickly", fontsize=8.5, color=INK)
axA.set_ylabel("action"); axA.grid(axis="y", color=GRID, lw=0.6)
axA.set_title("A   One trapped trajectory, n = 4, (T1,T4) = (2,8), γ = 0.1, start I2 = 200", loc="left", fontsize=10, color=INK)
plt.setp(axA.get_xticklabels(), visible=False)
axB.plot(t, np.cos(th3[s]), color=INK, lw=0.6)
axB.set_ylim(-1.15, 1.15); axB.set_yticks([-1, 0, 1]); axB.set_ylabel(r"$\cos\theta_3$")
axB.axvline(tb, color=MUTED, lw=0.8, ls=":"); axB.set_xlabel("t"); axB.grid(axis="y", color=GRID, lw=0.6)
axB.text(3, -0.55, r"$\cos\theta_3=-1$: $I_2$ and $I_3$ in anti-phase the whole time", fontsize=8.5, color=INK)
axB.set_title("B   Phase between the two middle modes (same trajectory)", loc="left", fontsize=10, color=INK)

# ---- C: action profile over the sites (mean over trajectories in the pair state at that time)
ramp = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
sites = np.arange(1, 5)
axC.plot([2], [200], color="#9a9a93", marker="o", ms=6, ls="none", label="t = 0: start, all in $I_2$")
for col, tt in zip(ramp, (1, 25, 50, 75, 100)):
    k = int(round(tt / 0.05)) - 1; sel = pair[:, k]
    axC.plot(sites, I[sel, k].mean(0), color=col, marker="o", ms=5, label=f"t = {tt}  ({sel.sum()} traj. still trapped)")
axC.plot(sites, prof_ness, color=INK, ls="--", marker="s", ms=4, lw=1.2, label="NESS (after break-up)")
axC.set_xticks(sites); axC.set_xticklabels([r"$I_1$" + "\n(T=2 bath)", r"$I_2$", r"$I_3$", r"$I_4$" + "\n(T=8 bath)"])
axC.set_ylabel("mean action"); axC.set_ylim(-3, 205); axC.grid(axis="y", color=GRID, lw=0.6)
axC.legend(fontsize=7.2, frameon=False, loc="upper right")
axC.set_title("C   Action profile while trapped (E0 = 200): two equal middle peaks that shrink", loc="left", fontsize=10, color=INK)

# ---- D: M(t), trapped trajectories, E0 = 100, 200, 400, vs n = 3
for col, E0 in zip(C[:3], (100, 200, 400)):
    a_, t_, I_, M_, th_, pair_, formed_, tau_ = trap_stats(E0)
    m = M_[formed_].mean(0); keep = t_ <= 330
    axD.plot(np.r_[0, t_[keep]], np.r_[E0, m[keep]], color=col)
    k5 = 99; ref = m[k5] - 0.78 * (t_ - t_[k5]); r = (t_ >= 5) & (ref > 0) & keep
    axD.plot(t_[r], ref[r], color=MUTED, ls=":", lw=1.1)
    axD.text(*{100: (40, 45), 200: (75, 115), 400: (165, 160)}[E0], f"E0={E0}", color=col, fontsize=9, fontweight="bold")
a3, t3 = load(f"{L}/drain/n3_E200_dump.f32", 3, 256, 0.002)
axD.plot(np.r_[0, t3], np.r_[200, a3[:, :, :3].sum(-1).mean(0)], color=C[3])
axD.text(120, 400, "blue/orange/green: n = 4, trapped trajectories\nyellow (near t = 0): n = 3, E0 = 200,\n    back to the NESS in ~2", color=INK, fontsize=8.3, va="top")
axD.text(120, 305, "dotted: slope −0.78 (constant leak)", color=MUTED, fontsize=8.3)
axD.set_xlim(0, 330); axD.set_ylim(0, 410); axD.set_xlabel("t"); axD.set_ylabel("mean M (trapped trajectories)")
axD.grid(axis="y", color=GRID, lw=0.6)
axD.set_title("D   Total action: straight-line decay, lifetime ∝ E0", loc="left", fontsize=10, color=INK)
fig.savefig(OUT, dpi=160, bbox_inches="tight", facecolor="white")
print("stream", s, "tau", tb)
