"""
LTE conjecture test for the resonant NLS energy cascade.

Two candidate forms are fit to log(rho) by linear regression:

  Form A:
    log rho ~ a + b * [ -(I_a^2 + I_b^2)/4 + I_a*I_b*cos(theta) ]

  Form B:
    log rho ~ a + b * [ (I_a^2 + I_b^2)/4 + I_a*I_b*(1 + cos(theta)) ]
    (Local Gibbs at temperature T_loc = -1/(2b))

  Form C (most general 2nd-order in I, with cos(theta) coupling):
    log rho ~ a + c1*I_a^2 + c2*I_b^2 + c3*I_a*I_b + c4*I_a*I_b*cos(theta)
    (Lets the data choose; check internal consistency)

Usage:
    python3 lte_analysis.py traj_N25.csv
"""

import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", help="Input CSV from lte_dump")
    ap.add_argument("--site", type=int, default=12,
                    help="Mode index j; we look at pair (j, j+1) (default 12)")
    ap.add_argument("--bins", type=int, default=40,
                    help="Histogram bins per axis (default 40)")
    ap.add_argument("--min-count", type=int, default=20,
                    help="Drop bins with fewer samples than this (default 20)")
    ap.add_argument("--I-max", type=float, default=None,
                    help="Cut I above this (default: 99th percentile)")
    args = ap.parse_args()

    print(f"Loading {args.csv} ...")
    df = pd.read_csv(args.csv)
    print(f"  {len(df)} samples, {df.shape[1]} columns")

    j = args.site
    if f"I_{j}" not in df.columns or f"I_{j+1}" not in df.columns:
        print(f"ERROR: columns I_{j} or I_{j+1} not in CSV")
        sys.exit(1)

    Ia = df[f"I_{j}"].to_numpy()
    Ib = df[f"I_{j+1}"].to_numpy()
    pa = df[f"phi_{j}"].to_numpy()
    pb = df[f"phi_{j+1}"].to_numpy()

    # theta_j = 2*(phi_{j+1} - phi_j), wrapped to [-pi, pi)
    theta = 2.0 * (pb - pa)
    theta = np.mod(theta + np.pi, 2*np.pi) - np.pi

    print(f"\nSite j={j} (pair {j},{j+1}):")
    print(f"  I_{j}: mean={Ia.mean():.3f}, std={Ia.std():.3f}, "
          f"min={Ia.min():.3e}, max={Ia.max():.3f}")
    print(f"  I_{j+1}: mean={Ib.mean():.3f}, std={Ib.std():.3f}, "
          f"min={Ib.min():.3e}, max={Ib.max():.3f}")
    print(f"  theta: mean={theta.mean():.3f}, "
          f"<cos>={np.cos(theta).mean():.4f}")

    # cap I to avoid ill-conditioned tails
    Imax = args.I_max if args.I_max else np.percentile(np.r_[Ia, Ib], 99)
    Imin = 1e-3
    print(f"\nUsing I in [{Imin:.3f}, {Imax:.3f}]")

    keep = (Ia > Imin) & (Ib > Imin) & (Ia < Imax) & (Ib < Imax)
    print(f"  kept {keep.sum()}/{len(Ia)} samples after I-cuts "
          f"({100*keep.mean():.1f}%)")

    Ia, Ib, theta = Ia[keep], Ib[keep], theta[keep]

    # --- 3D histogram in (I_a, I_b, theta) ---
    print(f"\nBuilding 3D histogram with {args.bins}^3 bins ...")
    Ia_edges = np.linspace(Imin, Imax, args.bins + 1)
    Ib_edges = np.linspace(Imin, Imax, args.bins + 1)
    th_edges = np.linspace(-np.pi, np.pi, args.bins + 1)

    H, edges = np.histogramdd(np.c_[Ia, Ib, theta],
                              bins=(Ia_edges, Ib_edges, th_edges))

    # bin centers
    Ia_c = 0.5*(Ia_edges[1:] + Ia_edges[:-1])
    Ib_c = 0.5*(Ib_edges[1:] + Ib_edges[:-1])
    th_c = 0.5*(th_edges[1:] + th_edges[:-1])

    # bin volumes (uniform in this case)
    dIa = Ia_edges[1] - Ia_edges[0]
    dIb = Ib_edges[1] - Ib_edges[0]
    dth = th_edges[1] - th_edges[0]
    vol = dIa * dIb * dth

    # density estimate
    rho = H / (H.sum() * vol)

    # mesh of bin centers
    IAm, IBm, THm = np.meshgrid(Ia_c, Ib_c, th_c, indexing='ij')
    counts = H

    # mask: bins with enough counts (drop low-density tails for fit)
    mask = counts >= args.min_count
    print(f"  bins above min-count={args.min_count}: "
          f"{mask.sum()}/{mask.size} ({100*mask.mean():.1f}%)")

    log_rho = np.log(rho[mask])
    Ia_m = IAm[mask]
    Ib_m = IBm[mask]
    th_m = THm[mask]

    # ============================================================
    # Form A:  log rho ~ a + b * F_A
    #          F_A = -(I_a^2 + I_b^2)/4 + I_a*I_b*cos(theta)
    # ============================================================
    F_A = -(Ia_m**2 + Ib_m**2)/4.0 + Ia_m*Ib_m*np.cos(th_m)
    A_design = np.c_[np.ones_like(F_A), F_A]
    coefs_A, res_A, _, _ = np.linalg.lstsq(A_design, log_rho, rcond=None)
    pred_A = A_design @ coefs_A
    ss_res_A = np.sum((log_rho - pred_A)**2)
    ss_tot = np.sum((log_rho - log_rho.mean())**2)
    R2_A = 1 - ss_res_A/ss_tot

    # ============================================================
    # Form B:  log rho ~ a + b * F_B
    #          F_B = (I_a^2 + I_b^2)/4 + I_a*I_b*(1 + cos(theta))
    # ============================================================
    F_B = (Ia_m**2 + Ib_m**2)/4.0 + Ia_m*Ib_m*(1.0 + np.cos(th_m))
    B_design = np.c_[np.ones_like(F_B), F_B]
    coefs_B, *_ = np.linalg.lstsq(B_design, log_rho, rcond=None)
    pred_B = B_design @ coefs_B
    ss_res_B = np.sum((log_rho - pred_B)**2)
    R2_B = 1 - ss_res_B/ss_tot

    # ============================================================
    # Form C: log rho ~ a + c1 I_a^2 + c2 I_b^2 + c3 I_a*I_b + c4 I_a*I_b*cos(theta)
    # ============================================================
    feats = np.c_[
        np.ones_like(Ia_m),
        Ia_m**2,
        Ib_m**2,
        Ia_m*Ib_m,
        Ia_m*Ib_m*np.cos(th_m),
    ]
    coefs_C, *_ = np.linalg.lstsq(feats, log_rho, rcond=None)
    pred_C = feats @ coefs_C
    ss_res_C = np.sum((log_rho - pred_C)**2)
    R2_C = 1 - ss_res_C/ss_tot

    print("\n" + "="*60)
    print(f"FIT RESULTS  (site j={j})")
    print("="*60)
    print(f"\nForm A:")
    print(f"  log rho = {coefs_A[0]:.4f} + {coefs_A[1]:.4f} * "
          f"[ -(I_a^2 + I_b^2)/4 + I_a I_b cos(theta) ]")
    print(f"  R^2 = {R2_A:.6f}")
    print(f"  Implied 'temperature' if read as -beta : "
          f"T = -1/(2*coefs_A[1]/some_scale) -- form needs interpretation")

    print(f"\nForm B:")
    print(f"  log rho = {coefs_B[0]:.4f} + {coefs_B[1]:.4f} * "
          f"[ (I_a^2 + I_b^2)/4 + I_a I_b (1 + cos(theta)) ]")
    print(f"  beta = -coefs_B[1] = {-coefs_B[1]:.4f}")
    if coefs_B[1] < 0:
        T_loc = -1.0 / (2.0 * coefs_B[1])
        print(f"  implied T_loc = -1/(2*coef) = {T_loc:.3f}  "
              f"(linear interp of T_1=10, T_n=2 at site j={j} -> "
              f"T = {10 - 8*j/24:.2f})")
    print(f"  R^2 = {R2_B:.6f}")

    print(f"\nForm C:")
    nm = ['const', 'I_a^2', 'I_b^2', 'I_a*I_b', 'I_a*I_b*cos']
    for n, c in zip(nm, coefs_C):
        print(f"  {n:18s}: {c:+.6f}")
    print(f"  R^2 = {R2_C:.6f}")
    print(f"\n  Internal consistency check (local Gibbs predicts):")
    print(f"    c1 == c2 (symmetry between I_a, I_b)?")
    print(f"      c1 = {coefs_C[1]:+.6f},  c2 = {coefs_C[2]:+.6f},  "
          f"|c1-c2|/|c1| = {abs(coefs_C[1]-coefs_C[2])/abs(coefs_C[1]):.3%}")
    print(f"    c3 == c4 (cross terms have same coefficient)?")
    print(f"      c3 = {coefs_C[3]:+.6f},  c4 = {coefs_C[4]:+.6f},  "
          f"ratio c3/c4 = {coefs_C[3]/coefs_C[4]:.4f}")
    print(f"    Form B predicts c1 = c2 = beta/4, c3 = c4 = beta")
    print(f"      => c1/c3 should equal 1/4 = 0.25")
    print(f"      observed c1/c3 = {coefs_C[1]/coefs_C[3]:.4f}")

    # ============================================================
    # Plot 1: predicted vs observed log rho
    # ============================================================
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, pred, R2, title in [
        (axes[0], pred_A, R2_A, "Form A"),
        (axes[1], pred_B, R2_B, "Form B"),
        (axes[2], pred_C, R2_C, "Form C"),
    ]:
        ax.scatter(pred, log_rho, s=2, alpha=0.3)
        lo = min(pred.min(), log_rho.min())
        hi = max(pred.max(), log_rho.max())
        ax.plot([lo, hi], [lo, hi], 'r--', lw=1)
        ax.set_xlabel("predicted log rho")
        ax.set_ylabel("observed log rho")
        ax.set_title(f"{title}\nR^2 = {R2:.5f}")
        ax.grid(alpha=0.3)
    plt.suptitle(f"LTE fit at site j={j} (pair {j},{j+1}), "
                 f"T_1=10, T_n=2, n=25", fontsize=12)
    plt.tight_layout()
    out1 = f"lte_fit_j{j}.png"
    plt.savefig(out1, dpi=120, bbox_inches='tight')
    print(f"\nSaved: {out1}")

    # ============================================================
    # Plot 2: 1D marginals of (I_a, I_b, theta) for visual sanity
    # ============================================================
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    axes[0].hist(Ia, bins=80, density=True, alpha=0.6, color='C0')
    axes[0].set_xlabel(f"I_{j}"); axes[0].set_ylabel("p")
    axes[0].set_title(f"marginal I_{j}")
    axes[1].hist(Ib, bins=80, density=True, alpha=0.6, color='C1')
    axes[1].set_xlabel(f"I_{j+1}")
    axes[1].set_title(f"marginal I_{j+1}")
    axes[2].hist(theta, bins=80, density=True, alpha=0.6, color='C2')
    axes[2].set_xlabel(f"theta_{j} = 2(phi_{j+1}-phi_{j})")
    axes[2].set_title("marginal theta")
    axes[2].axvline(0, color='k', lw=0.5)
    plt.suptitle(f"1D marginals at site j={j}", fontsize=12)
    plt.tight_layout()
    out2 = f"lte_marginals_j{j}.png"
    plt.savefig(out2, dpi=120, bbox_inches='tight')
    print(f"Saved: {out2}")

    # ============================================================
    # Plot 3: 2D slices  log rho(I_a, I_b) at theta=0,
    #                    log rho(I_a, theta) at I_b ~ <I_b>
    # ============================================================
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # slice in theta=0 (the central theta bin)
    th_idx = args.bins // 2
    slab = rho[:, :, th_idx]
    with np.errstate(divide='ignore'):
        logslab = np.log(slab + 1e-30)
    im0 = axes[0].pcolormesh(Ia_c, Ib_c, logslab.T,
                             cmap='viridis', shading='auto')
    axes[0].set_xlabel(f"I_{j}"); axes[0].set_ylabel(f"I_{j+1}")
    axes[0].set_title(f"log rho at theta ~ 0")
    plt.colorbar(im0, ax=axes[0])

    # slice in I_b ~ <I_b>
    Ib_idx = np.argmin(np.abs(Ib_c - Ib.mean()))
    slab2 = rho[:, Ib_idx, :]
    with np.errstate(divide='ignore'):
        logslab2 = np.log(slab2 + 1e-30)
    im1 = axes[1].pcolormesh(Ia_c, th_c, logslab2.T,
                             cmap='viridis', shading='auto')
    axes[1].set_xlabel(f"I_{j}"); axes[1].set_ylabel("theta")
    axes[1].set_title(f"log rho at I_{j+1} ~ {Ib_c[Ib_idx]:.2f}")
    plt.colorbar(im1, ax=axes[1])

    plt.suptitle(f"2D slices of joint density, site j={j}", fontsize=12)
    plt.tight_layout()
    out3 = f"lte_slices_j{j}.png"
    plt.savefig(out3, dpi=120, bbox_inches='tight')
    print(f"Saved: {out3}")

    print("\nDone.")

if __name__ == "__main__":
    main()
