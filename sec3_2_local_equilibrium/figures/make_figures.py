"""Fig. 5 of Sec. 3.2: residual of the NESS histogram against the rescaled equilibrium histogram, hot side.

Usage: python3 make_figures.py [OUT]      (OUT: folder for the PDF; default: this folder)
Reads ../results/aggregated_hist_v2.npz (80^3 histograms of (I_a, I_{a+1}, theta_{a+1}), summed over all batches,
for the NESS (10,5) and the equilibrium T = 7.5 runs) and ../results/lte_summary_v2.csv (the fitted slopes), and
writes lte_residual_mesh_hot_side.pdf.
For each bond: on the bins with I_a, I_{a+1} < 2.5 and at least 50 counts in both histograms, fit
log q = b0 + b1 log p (q = NESS, p = equilibrium, weights q + p); the residual log q - (b0 + b1 log p) is averaged
over theta with weights q + p and drawn as a mesh over (I_a, I_{a+1}) (MATLAB mesh() style: white faces, edges
coloured by height).  The slope b1 is checked against the analysis table.
"""
import sys, os, csv, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "..", "common"))
import paperstyle as S
import matplotlib as mpl, matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
AGG = os.path.join(HERE, "..", "results", "aggregated_hist_v2.npz"); SUM = os.path.join(HERE, "..", "results", "lte_summary_v2.csv")
OUT = os.path.join(sys.argv[1] if len(sys.argv) > 1 else HERE, "lte_residual_mesh_hot_side.pdf")
NB, CORE, MINC = 80, 2.5, 50
centers = np.linspace(0.025, 3.975, NB)
agg = np.load(AGG); summ = {(int(r["n"]), int(r["site_left_0based"])): r for r in csv.DictReader(open(SUM))}
PARULA = LinearSegmentedColormap.from_list("parula_like", ["#352a87", "#0f5cdd", "#1481d6", "#06a4ca", "#2eb7a4", "#87bf77", "#d1bb59", "#fec832", "#f9fb0e"])
def residual_surface(q, p):
    q = q.astype(float); p = p.astype(float)
    core2 = (centers[:, None] < CORE) & (centers[None, :] < CORE)
    core = np.broadcast_to(core2[:, :, None], q.shape)
    reliable = (q >= MINC) & (p >= MINC) & core
    qc, pc = q[core], p[core]; rel = (qc >= MINC) & (pc >= MINC)
    x = np.log(pc[rel] / pc.sum()); y = np.log(qc[rel] / qc.sum()); w = qc[rel] + pc[rel]
    A = np.column_stack([np.ones_like(x), x]) * np.sqrt(w)[:, None]
    b0, b1 = np.linalg.lstsq(A, y * np.sqrt(w), rcond=None)[0]
    res = np.zeros_like(q)
    res[reliable] = np.log(q[reliable] / q.sum()) - (b0 + b1 * np.log(p[reliable] / p.sum()))
    pw = (q + p) * reliable; num = (res * pw).sum(2); den = pw.sum(2)
    return np.where(den > 0, num / np.where(den > 0, den, 1), np.nan), b1
keys = [(25, 6), (50, 12), (100, 24)]; surf = {}
for n, site in keys:
    surf[(n, site)], b1 = residual_surface(agg[f"NESS_n{n}_j{site}"], agg[f"EQ_n{n}_j{site}"])
    assert abs(b1 - float(summ[(n, site)]["slope_x"])) < 1e-6
scale = np.percentile(np.concatenate([np.abs(s[np.isfinite(s)]) for s in surf.values()]), 98)
zall = np.concatenate([s[np.ix_(centers < CORE, centers < CORE)][np.isfinite(s[np.ix_(centers < CORE, centers < CORE)])] for s in surf.values()])
zlo, zhi = zall.min(), zall.max(); print('zrange', zlo, zhi)
fig = plt.figure(figsize=(S.TEXTWIDTH, 2.6)); L = S.letters()
m = centers < CORE; X, Y = np.meshgrid(centers[m], centers[m])
for i, (n, site) in enumerate(keys):
    ax = fig.add_subplot(1, 3, i + 1, projection="3d")
    Z = surf[(n, site)].T[np.ix_(m, m)]
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    polys, vals = [], []
    for r in range(Z.shape[0] - 1):
        for c in range(Z.shape[1] - 1):
            q = [(r, c), (r, c + 1), (r + 1, c + 1), (r + 1, c)]
            zz = [Z[i, j] for i, j in q]
            if np.all(np.isfinite(zz)):
                polys.append([(X[i, j], Y[i, j], Z[i, j]) for i, j in q]); vals.append(np.mean(zz))
    pc = Poly3DCollection(polys, facecolors="w", edgecolors=PARULA(np.clip(mpl.colors.Normalize(-scale, scale)(np.array(vals)), 0, 1)), linewidths=0.3)
    ax.add_collection3d(pc)
    ax.set_xlim(0, CORE); ax.set_ylim(0, CORE); ax.set_zlim(zlo, zhi)
    ax.view_init(elev=30, azim=-127.5); ax.set_xticks([0.5, 1, 1.5, 2, 2.5]); ax.set_yticks([0.5, 1, 1.5, 2, 2.5])
    ax.set_xlabel(r"$I_a$", labelpad=-8); ax.set_ylabel(r"$I_{a+1}$", labelpad=-8)
    ax.tick_params(axis="both", which="major", pad=-3, labelsize=7)
    for p in (ax.xaxis, ax.yaxis, ax.zaxis):
        p.pane.fill = False; p.pane.set_edgecolor("0.75"); p._axinfo["grid"].update(color="0.88", linewidth=0.4)
    ax.set_title(f"({L[i]}) $n={n}$", y=1.0, pad=0)
fig.subplots_adjust(left=0.0, right=0.98, wspace=0.02, bottom=0.05, top=0.95)
fig.savefig(OUT)
print("scale", scale)
