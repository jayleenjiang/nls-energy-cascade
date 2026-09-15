"""
Density mesh visualizer for .hist files produced by lte_histogram.cpp.

Produces several visualizations of rho(I_a, I_b, theta) from one .hist file:
  1. 1D marginals: P(I_a), P(I_b), P(theta)
  2. 2D marginals: rho(I_a, I_b), rho(I_a, theta), rho(I_b, theta)
  3. 2D theta slices: rho(I_a, I_b | theta=fixed) for 4 theta values
  4. 3D scatter mesh of high-density voxels

Usage:
    python3 hist_visualize.py histo_N25_j12.hist
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
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
    ap.add_argument("--I-max-display", type=float, default=3.0,
                    help="cap I axis display range (default 3.0)")
    args = ap.parse_args()

    H = load_hist(args.histfile)
    NB = H["NB"]
    counts = H["counts"]
    I_LO, I_HI = H["I_LO"], H["I_HI"]
    TH_LO, TH_HI = H["TH_LO"], H["TH_HI"]
    total = counts.sum()

    print(f"Loaded {args.histfile}: NB={NB}, total samples = {total:,}")

    dI = (I_HI - I_LO) / NB
    dth = (TH_HI - TH_LO) / NB
    Ia_c = I_LO + dI * (np.arange(NB) + 0.5)
    Ib_c = I_LO + dI * (np.arange(NB) + 0.5)
    th_c = TH_LO + dth * (np.arange(NB) + 0.5)

    bin_vol = dI * dI * dth
    rho_3d = counts / (total * bin_vol)

    # ============================================================
    # 1. 1D marginals
    # ============================================================
    P_Ia = rho_3d.sum(axis=(1, 2)) * dI * dth  # integrate over b, theta
    P_Ib = rho_3d.sum(axis=(0, 2)) * dI * dth
    P_th = rho_3d.sum(axis=(0, 1)) * dI * dI

    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].plot(Ia_c, P_Ia, lw=2)
    axes[0].set_xlabel(r"$I_{12}$")
    axes[0].set_ylabel(r"$P(I_{12})$")
    axes[0].set_title(r"1D marginal $P(I_{12})$")
    axes[0].set_xlim(0, args.I_max_display)
    axes[0].grid(alpha=0.3)

    axes[1].plot(Ib_c, P_Ib, lw=2, color='C1')
    axes[1].set_xlabel(r"$I_{13}$")
    axes[1].set_ylabel(r"$P(I_{13})$")
    axes[1].set_title(r"1D marginal $P(I_{13})$")
    axes[1].set_xlim(0, args.I_max_display)
    axes[1].grid(alpha=0.3)

    axes[2].plot(th_c, P_th, lw=2, color='C2')
    axes[2].set_xlabel(r"$\theta_{12}$")
    axes[2].set_ylabel(r"$P(\theta_{12})$")
    axes[2].set_title(r"1D marginal $P(\theta_{12})$")
    axes[2].axvline(0, color='k', lw=0.5)
    axes[2].grid(alpha=0.3)

    plt.suptitle(f"1D marginals ({total:,} samples)", fontsize=13)
    plt.tight_layout()
    out1 = args.histfile.replace(".hist", "_1d.png")
    plt.savefig(out1, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out1}")

    # ============================================================
    # 2. 2D marginals
    # ============================================================
    rho_IaIb = rho_3d.sum(axis=2) * dth          # integrate over theta
    rho_Iath = rho_3d.sum(axis=1) * dI           # integrate over b
    rho_Ibth = rho_3d.sum(axis=0) * dI           # integrate over a

    fig, axes = plt.subplots(1, 3, figsize=(16, 5))

    # rho(I_a, I_b)
    rho_show = np.where(rho_IaIb > 0, rho_IaIb, np.nan)
    im0 = axes[0].pcolormesh(Ia_c, Ib_c, rho_show.T,
                             cmap='viridis', shading='auto',
                             norm=LogNorm(
                                 vmin=max(np.nanmin(rho_show), 1e-4),
                                 vmax=np.nanmax(rho_show)))
    axes[0].set_xlabel(r"$I_{12}$")
    axes[0].set_ylabel(r"$I_{13}$")
    axes[0].set_title(r"2D marginal $\rho(I_{12}, I_{13})$")
    axes[0].set_aspect('equal')
    axes[0].set_xlim(0, args.I_max_display)
    axes[0].set_ylim(0, args.I_max_display)
    plt.colorbar(im0, ax=axes[0], label=r'$\rho$')

    # rho(I_a, theta)
    rho_show = np.where(rho_Iath > 0, rho_Iath, np.nan)
    im1 = axes[1].pcolormesh(Ia_c, th_c, rho_show.T,
                             cmap='viridis', shading='auto',
                             norm=LogNorm(
                                 vmin=max(np.nanmin(rho_show), 1e-4),
                                 vmax=np.nanmax(rho_show)))
    axes[1].set_xlabel(r"$I_{12}$")
    axes[1].set_ylabel(r"$\theta_{12}$")
    axes[1].set_title(r"2D marginal $\rho(I_{12}, \theta_{12})$")
    axes[1].set_xlim(0, args.I_max_display)
    plt.colorbar(im1, ax=axes[1], label=r'$\rho$')

    # rho(I_b, theta)
    rho_show = np.where(rho_Ibth > 0, rho_Ibth, np.nan)
    im2 = axes[2].pcolormesh(Ib_c, th_c, rho_show.T,
                             cmap='viridis', shading='auto',
                             norm=LogNorm(
                                 vmin=max(np.nanmin(rho_show), 1e-4),
                                 vmax=np.nanmax(rho_show)))
    axes[2].set_xlabel(r"$I_{13}$")
    axes[2].set_ylabel(r"$\theta_{12}$")
    axes[2].set_title(r"2D marginal $\rho(I_{13}, \theta_{12})$")
    axes[2].set_xlim(0, args.I_max_display)
    plt.colorbar(im2, ax=axes[2], label=r'$\rho$')

    plt.suptitle(f"2D marginals (one variable integrated out)",
                 fontsize=13)
    plt.tight_layout()
    out2 = args.histfile.replace(".hist", "_2d.png")
    plt.savefig(out2, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out2}")

    # ============================================================
    # 3. 2D slices at fixed theta
    # ============================================================
    theta_targets = [-np.pi, -np.pi/2, 0.0, np.pi/2]
    labels = [r"$\theta = -\pi$", r"$\theta = -\pi/2$",
              r"$\theta = 0$", r"$\theta = +\pi/2$"]
    fig, axes = plt.subplots(1, 4, figsize=(20, 5))
    for ax, theta0, label in zip(axes, theta_targets, labels):
        # nearest bin in theta
        it = np.argmin(np.abs(th_c - theta0))
        slab = counts[:, :, it].astype(float)
        # normalize to a conditional density (within the slab)
        s = slab.sum()
        if s > 0:
            slab_dens = slab / (s * dI * dI)
        else:
            slab_dens = slab

        slab_show = np.where(slab_dens > 0, slab_dens, np.nan)
        if not np.all(np.isnan(slab_show)):
            vmin = max(np.nanmin(slab_show), 1e-4)
            vmax = np.nanmax(slab_show)
            im = ax.pcolormesh(Ia_c, Ib_c, slab_show.T,
                               cmap='viridis', shading='auto',
                               norm=LogNorm(vmin=vmin, vmax=vmax))
            plt.colorbar(im, ax=ax, label=r'$\rho(\cdot|\theta)$')

        ax.set_xlabel(r"$I_{12}$")
        ax.set_ylabel(r"$I_{13}$")
        ax.set_title(label + f"\n({int(slab.sum()):,} samples in slab)")
        ax.set_aspect('equal')
        ax.set_xlim(0, args.I_max_display)
        ax.set_ylim(0, args.I_max_display)

    plt.suptitle(r"Conditional density $\rho(I_{12}, I_{13}\,|\,\theta)$ "
                 r"at fixed $\theta$", fontsize=13)
    plt.tight_layout()
    out3 = args.histfile.replace(".hist", "_theta_slices.png")
    plt.savefig(out3, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out3}")

    # ============================================================
    # 4. 3D scatter of high-density voxels
    # ============================================================
    # pick bins with density above threshold (top X% by count)
    thresh_count = np.percentile(counts[counts > 0], 90)  # top 10%
    mask = counts >= thresh_count
    ia, ib, it = np.where(mask)
    a_v = Ia_c[ia]
    b_v = Ib_c[ib]
    t_v = th_c[it]
    c_v = counts[ia, ib, it].astype(float)

    fig = plt.figure(figsize=(14, 6))

    # left: log color, full theta range
    ax = fig.add_subplot(1, 2, 1, projection='3d')
    p = ax.scatter(a_v, b_v, t_v, c=np.log(c_v),
                   cmap='viridis', s=8, alpha=0.5)
    ax.set_xlabel(r"$I_{12}$")
    ax.set_ylabel(r"$I_{13}$")
    ax.set_zlabel(r"$\theta_{12}$")
    ax.set_title(f"High-density voxels (top 10%)\n"
                 f"color = log(count), {mask.sum():,} voxels")
    ax.set_xlim(0, args.I_max_display)
    ax.set_ylim(0, args.I_max_display)
    plt.colorbar(p, ax=ax, shrink=0.6, label='log(count)')

    # right: side view, peak structure
    ax2 = fig.add_subplot(1, 2, 2, projection='3d')
    # show only the very densest bins
    thresh_top = np.percentile(counts[counts > 0], 99)
    mask_top = counts >= thresh_top
    ia2, ib2, it2 = np.where(mask_top)
    p2 = ax2.scatter(Ia_c[ia2], Ib_c[ib2], th_c[it2],
                     c=counts[ia2, ib2, it2],
                     cmap='hot', s=15, alpha=0.7)
    ax2.set_xlabel(r"$I_{12}$")
    ax2.set_ylabel(r"$I_{13}$")
    ax2.set_zlabel(r"$\theta_{12}$")
    ax2.set_title(f"Peak voxels (top 1%)\n{mask_top.sum():,} voxels")
    ax2.set_xlim(0, args.I_max_display)
    ax2.set_ylim(0, args.I_max_display)
    plt.colorbar(p2, ax=ax2, shrink=0.6, label='count')

    plt.suptitle(f"3D mesh of high-density voxels", fontsize=13)
    plt.tight_layout()
    out4 = args.histfile.replace(".hist", "_3d.png")
    plt.savefig(out4, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out4}")

    print("\nDone. Four files:")
    print(f"  {out1}  -- 1D marginals P(I_a), P(I_b), P(theta)")
    print(f"  {out2}  -- 2D marginals integrated over one axis")
    print(f"  {out3}  -- 2D slices at fixed theta")
    print(f"  {out4}  -- 3D voxel mesh")


if __name__ == "__main__":
    main()
