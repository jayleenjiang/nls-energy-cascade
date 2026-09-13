"""
Fixed-theta 3D mesh visualizer.

Per advisor: pick several theta values, slice the histogram at fixed
theta, plot -log P(I_1, I_2 | theta) as a 3D surface mesh.  In
particular, at cos(theta) = 0 (i.e. theta = +-pi/2) a strict Gibbs
density would give -log P = const + (I_1^2 + I_2^2), i.e. a circular
paraboloid with no cross term.

Usage:
    python3 theta_slice_mesh.py cond_S12.5_j12.hist

Output:
    one PNG per theta value containing a 3D mesh of -log P over
    (I_1, I_2).
"""

import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def load_hist(fname):
    meta = {}
    with open(fname) as f:
        lines = f.readlines()
    idx = 0
    for i, ln in enumerate(lines):
        ln = ln.strip()
        if ln.startswith("#") or ln == "":
            continue
        parts = ln.split()
        if parts[0] == "NB":
            meta["NB"] = int(parts[1])
        elif parts[0] in ("I_LO", "I_HI", "TH_LO", "TH_HI"):
            meta[parts[0]] = float(parts[1])
        elif parts[0] in ("TOTAL", "OVERFLOW"):
            meta[parts[0]] = int(parts[1])
        elif len(parts) == 4:
            idx = i
            break
    NB = meta["NB"]
    counts = np.zeros((NB, NB, NB), dtype=np.int64)
    for ln in lines[idx:]:
        ln = ln.strip()
        if not ln or ln.startswith("#"):
            continue
        p = ln.split()
        if len(p) != 4:
            continue
        ia, ib, it, c = int(p[0]), int(p[1]), int(p[2]), int(p[3])
        counts[ia, ib, it] = c
    meta["counts"] = counts
    return meta


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("histfile")
    ap.add_argument("--I-max", type=float, default=2.5,
                    help="cap I axis range for plot (default 2.5)")
    ap.add_argument("--min-count", type=int, default=20,
                    help="bins with fewer counts shown as NaN (default 20)")
    ap.add_argument("--n-slices", type=int, default=5,
                    help="number of theta slices to plot (default 5)")
    args = ap.parse_args()

    H = load_hist(args.histfile)
    NB = H["NB"]
    counts = H["counts"]
    I_LO, I_HI = H["I_LO"], H["I_HI"]
    TH_LO, TH_HI = H["TH_LO"], H["TH_HI"]
    total = H.get("TOTAL", counts.sum())
    print(f"Loaded {args.histfile}: NB={NB}, total={total:,}")

    dI = (I_HI - I_LO) / NB
    dth = (TH_HI - TH_LO) / NB
    Ia_c = I_LO + dI * (np.arange(NB) + 0.5)
    Ib_c = I_LO + dI * (np.arange(NB) + 0.5)
    th_c = TH_LO + dth * (np.arange(NB) + 0.5)

    # theta values to slice at
    # include cos(theta) = 0 explicitly (i.e. theta = pi/2)
    target_thetas = [-np.pi + dth, -np.pi/2, 0.0, np.pi/2, np.pi - dth][:args.n_slices]
    th_idx = [int(np.argmin(np.abs(th_c - t))) for t in target_thetas]

    # restrict to I < I_max for cleaner plots
    Ia_keep = Ia_c < args.I_max
    Ib_keep = Ib_c < args.I_max
    Ia_show = Ia_c[Ia_keep]
    Ib_show = Ib_c[Ib_keep]
    nA, nB = len(Ia_show), len(Ib_show)

    Iam, Ibm = np.meshgrid(Ia_show, Ib_show, indexing='ij')

    # for each theta slab compute -log P_slab(I_a, I_b)
    # we use a single normalization per slab
    fig = plt.figure(figsize=(5 * len(th_idx), 5))

    base = os.path.basename(args.histfile).replace(".hist", "")

    for k, (it, theta_val) in enumerate(zip(th_idx, target_thetas)):
        slab = counts[:, :, it].astype(float)         # NB x NB
        # restrict to display range
        slab = slab[np.ix_(Ia_keep, Ib_keep)]         # nA x nB
        # mask sparse bins
        keep = slab >= args.min_count
        # normalize within slab (conditional pdf in I_a, I_b for fixed slab)
        slab_total = slab.sum()
        if slab_total > 0:
            p = slab / (slab_total * dI * dI)         # conditional density
        else:
            p = slab
        neg_logp = np.full_like(p, np.nan, dtype=float)
        neg_logp[keep] = -np.log(p[keep])

        ax = fig.add_subplot(1, len(th_idx), k + 1, projection='3d')
        # surface: x = Iam, y = Ibm, z = neg_logp
        # use np.ma so NaN doesn't break plot_surface
        nlp_masked = np.where(np.isnan(neg_logp), np.nanmax(neg_logp) + 0.5, neg_logp)
        surf = ax.plot_surface(Iam, Ibm, nlp_masked,
                               cmap='viridis',
                               linewidth=0.2, antialiased=True,
                               rstride=2, cstride=2, edgecolor='none')
        ax.set_xlabel(r"$I_a$")
        ax.set_ylabel(r"$I_b$")
        ax.set_zlabel(r"$-\log P$")
        cos_val = np.cos(theta_val)
        ax.set_title(
            rf"$\theta={theta_val:+.2f}$,  $\cos\theta={cos_val:+.2f}$" + "\n"
            + f"({int(slab_total):,} samples in slab)",
            fontsize=10)

    plt.suptitle(f"{base}:  $-\\log P(I_a, I_b\\,|\\,\\theta)$ at fixed $\\theta$",
                 fontsize=12)
    plt.tight_layout()
    out = args.histfile.replace(".hist", "_thetamesh.png")
    plt.savefig(out, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved 3D mesh: {out}")

    # ---------- supplementary: at cos theta = 0, overlay Gibbs paraboloid ----
    # find the slice closest to theta = pi/2 (cos=0)
    pi2_idx = int(np.argmin(np.abs(th_c - np.pi/2)))
    slab = counts[:, :, pi2_idx].astype(float)
    slab = slab[np.ix_(Ia_keep, Ib_keep)]
    keep = slab >= args.min_count
    slab_total = slab.sum()
    if slab_total > 0:
        p = slab / (slab_total * dI * dI)
    else:
        p = slab
    neg_logp = np.full_like(p, np.nan, dtype=float)
    neg_logp[keep] = -np.log(p[keep])

    # naive Gibbs prediction at cos theta = 0:
    #   strict-Gibbs (Form B):     -log P propto (I_a^2 + I_b^2)/4 + I_a I_b
    #   conjecture (Form A) at 0:  -log P propto (I_a^2 + I_b^2)/4   (no cross term)
    # We overlay both, anchored at the same constant offset (chosen so the
    # surfaces match the data at the minimum of the data).
    data_min = np.nanmin(neg_logp)
    formA_pred = (Iam**2 + Ibm**2) / 4.0
    formA_pred = formA_pred - formA_pred.min() + data_min
    formB_pred = (Iam**2 + Ibm**2) / 4.0 + Iam * Ibm * (1 + np.cos(np.pi/2))
    formB_pred = formB_pred - formB_pred.min() + data_min

    fig = plt.figure(figsize=(15, 5))
    titles = [
        f"data: $-\\log P(I_a, I_b\\,|\\,\\theta=\\pi/2)$",
        f"Form A (conjecture) at $\\cos\\theta=0$:\n"
        r"$\propto (I_a^2 + I_b^2)/4$",
        f"Form B (strict local Gibbs) at $\\cos\\theta=0$:\n"
        r"$\propto (I_a^2 + I_b^2)/4 + I_a I_b$",
    ]
    surfaces = [neg_logp, formA_pred, formB_pred]
    for k in range(3):
        ax = fig.add_subplot(1, 3, k + 1, projection='3d')
        Z = surfaces[k]
        Z_plot = np.where(np.isnan(Z), np.nanmax(Z) + 0.5, Z) if k == 0 else Z
        ax.plot_surface(Iam, Ibm, Z_plot, cmap='viridis',
                        linewidth=0.2, antialiased=True,
                        rstride=2, cstride=2, edgecolor='none')
        ax.set_xlabel(r"$I_a$")
        ax.set_ylabel(r"$I_b$")
        ax.set_zlabel(r"$-\log P$")
        ax.set_title(titles[k], fontsize=10)

    plt.suptitle(f"{base}:  data vs Form A vs Form B at $\\cos\\theta=0$",
                 fontsize=12)
    plt.tight_layout()
    out2 = args.histfile.replace(".hist", "_costheta0_compare.png")
    plt.savefig(out2, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved cos(theta)=0 comparison: {out2}")
    print("Done.")


if __name__ == "__main__":
    main()
