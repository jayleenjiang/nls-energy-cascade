"""
Form A:  log rho ~ -1/4(a^2 + b^2) + a*b*cos(theta)
Form B:  log rho ~ -1/4(a^2 + b^2) - a*b*(1 + cos(theta))

Usage:
    python3 theta_slice_analysis.py traj_N25.csv
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.lines import Line2D


def make_empirical_figure(theta_slices, rho_panels, I_c, j, delta, out):
    n_cols = len(theta_slices)
    fig, axes = plt.subplots(1, n_cols, figsize=(5*n_cols, 5))
    if n_cols == 1:
        axes = [axes]
    for col, (label, theta0, rho, log_rho, n_in) in enumerate(rho_panels):
        ax = axes[col]
        im = ax.pcolormesh(I_c, I_c, log_rho,
                           cmap='viridis', shading='auto',
                           vmin=np.log(1e-3), vmax=log_rho.max())
        ax.set_xlabel(f"I_{j}")
        ax.set_ylabel(f"I_{j+1}")
        ax.set_title(f"{label}\nempirical log rho ({n_in:,} samples)")
        ax.set_aspect('equal')
        plt.colorbar(im, ax=ax)
    plt.suptitle(
        f"Empirical conditional density rho(I_{j}, I_{j+1} | theta)\n"
        f"theta half-width = {delta:.2f} rad", fontsize=13)
    plt.tight_layout()
    plt.savefig(out, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def make_model_figure(theta_slices, rho_panels, I_c, j, delta,
                      Imax, A, B, model_name, color, ls, out):
    n_cols = len(theta_slices)
    fig, axes = plt.subplots(1, n_cols, figsize=(5*n_cols, 5))
    if n_cols == 1:
        axes = [axes]
    for col, (label, theta0, rho, log_rho, n_in) in enumerate(rho_panels):
        ax = axes[col]

        rho_show = np.where(rho.T > 0, rho.T, np.nan)
        if np.nanmax(rho_show) > 0:
            im = ax.pcolormesh(
                I_c, I_c, rho_show,
                cmap='Greys',
                norm=LogNorm(
                    vmin=max(rho_show[rho_show > 0].min(), 1e-4),
                    vmax=np.nanmax(rho_show)),
                shading='auto')
            plt.colorbar(im, ax=ax, label='rho (empirical)')

        cos_t = np.cos(theta0)
        if model_name == "Form A":
            F_model = -0.25*(A**2 + B**2) + A*B*cos_t
        else:  # Form B
            F_model = -0.25*(A**2 + B**2) - A*B*(1.0 + cos_t)

        if np.isfinite(log_rho).any():
            lo = np.nanpercentile(log_rho[log_rho > -np.inf], 50)
            hi = np.nanpercentile(log_rho[log_rho > -np.inf], 99)
        else:
            lo, hi = -10, 0

        levels = np.linspace(F_model.max() - (hi - lo), F_model.max(), 6)
        ax.contour(A, B, F_model, levels=levels,
                   colors=color, linewidths=1.4, linestyles=ls)

        legend_elems = [Line2D([0], [0], color=color, ls=ls,
                               label=f"{model_name} contours")]
        ax.legend(handles=legend_elems, loc='upper right', fontsize=9)

        ax.set_xlabel(f"I_{j}")
        ax.set_ylabel(f"I_{j+1}")
        ax.set_title(f"{label}\nempirical (gray) + {model_name}")
        ax.set_aspect('equal')
        ax.set_xlim(0, Imax)
        ax.set_ylim(0, Imax)

    plt.suptitle(
        f"{model_name} contours overlaid on empirical density\n"
        f"theta half-width = {delta:.2f} rad", fontsize=13)
    plt.tight_layout()
    plt.savefig(out, dpi=120, bbox_inches='tight')
    plt.close(fig)
    print(f"Saved: {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--site", type=int, default=12,
                    help="mode index j; analyze pair (j, j+1)")
    ap.add_argument("--bins", type=int, default=60)
    ap.add_argument("--I-max", type=float, default=None,
                    help="cap I (default: 99th percentile)")
    ap.add_argument("--theta-half-width", type=float, default=0.2,
                    help="half-width of theta slice in radians (default 0.2)")
    args = ap.parse_args()

    j = args.site
    print(f"Loading {args.csv} ...")
    df = pd.read_csv(args.csv)
    print(f"  {len(df)} samples")

    Ia = df[f"I_{j}"].to_numpy()
    Ib = df[f"I_{j+1}"].to_numpy()
    pa = df[f"phi_{j}"].to_numpy()
    pb = df[f"phi_{j+1}"].to_numpy()
    theta = 2.0 * (pb - pa)
    theta = np.mod(theta + np.pi, 2*np.pi) - np.pi

    Imax = args.I_max if args.I_max else np.percentile(np.r_[Ia, Ib], 99)
    Imin = 0.0
    print(f"Using I in [{Imin:.3f}, {Imax:.3f}]")

    theta_slices = [
        ("theta = -pi",   -np.pi),
        ("theta = -pi/2", -np.pi/2),
        ("theta = 0",      0.0),
        ("theta = +pi/2",  np.pi/2),
    ]

    delta = args.theta_half_width

    Igrid = np.linspace(0.01, Imax, 200)
    A, B = np.meshgrid(Igrid, Igrid)

    I_edges = np.linspace(Imin, Imax, args.bins + 1)
    I_c = 0.5*(I_edges[1:] + I_edges[:-1])
    bin_area = (I_edges[1] - I_edges[0])**2

    print(f"\n{'theta':>12s}  {'mask wrap':>12s}  {'samples':>10s}  "
          f"{'<a>':>8s}  {'<b>':>8s}")

    rho_panels = []
    for label, theta0 in theta_slices:
        diff = theta - theta0
        diff = np.mod(diff + np.pi, 2*np.pi) - np.pi
        mask = np.abs(diff) <= delta
        n_in = mask.sum()
        a_in, b_in = Ia[mask], Ib[mask]

        print(f"{label:>12s}  +/- {delta:.2f}  {n_in:>10d}  "
              f"{a_in.mean():8.3f}  {b_in.mean():8.3f}")

        H, _, _ = np.histogram2d(a_in, b_in, bins=[I_edges, I_edges])
        rho = H / (n_in * bin_area + 1e-12)

        with np.errstate(divide='ignore'):
            log_rho = np.log(rho.T + 1e-12)

        rho_panels.append((label, theta0, rho, log_rho, n_in))

    # ============================================================
    # Three separate figure files
    # ============================================================
    make_empirical_figure(
        theta_slices, rho_panels, I_c, j, delta,
        out=f"theta_slices_empirical_j{j}.png")

    make_model_figure(
        theta_slices, rho_panels, I_c, j, delta,
        Imax, A, B,
        model_name="Form A", color="red",  ls="--",
        out=f"theta_slices_formA_j{j}.png")

    make_model_figure(
        theta_slices, rho_panels, I_c, j, delta,
        Imax, A, B,
        model_name="Form B", color="blue", ls="-",
        out=f"theta_slices_formB_j{j}.png")

    # ============================================================
    # Quantitative fit
    # ============================================================
    print("\n=== Fitting log rho = c0 + c1*(a^2+b^2) + c2*ab on each theta slice ===")
    print(f"{'theta':>12s}  {'c0':>8s}  {'c1':>8s}  {'c2':>8s}  "
          f"{'A pred c2':>10s}  {'B pred c2':>10s}  {'R^2':>6s}")

    for label, theta0 in theta_slices:
        diff = theta - theta0
        diff = np.mod(diff + np.pi, 2*np.pi) - np.pi
        mask = np.abs(diff) <= delta
        a_in, b_in = Ia[mask], Ib[mask]

        keep = (a_in > 0.01) & (b_in > 0.01) & (a_in < Imax) & (b_in < Imax)
        a_in, b_in = a_in[keep], b_in[keep]

        H, _, _ = np.histogram2d(a_in, b_in, bins=[I_edges, I_edges])
        AA, BB = np.meshgrid(I_c, I_c, indexing='ij')
        dens = H / (H.sum() * bin_area + 1e-12)
        m = H >= 20
        if m.sum() < 20:
            print(f"{label:>12s}  insufficient bins")
            continue
        log_rho_m = np.log(dens[m])
        a_m = AA[m]
        b_m = BB[m]
        feats = np.c_[np.ones_like(a_m), a_m**2 + b_m**2, a_m*b_m]
        coefs, *_ = np.linalg.lstsq(feats, log_rho_m, rcond=None)
        pred = feats @ coefs
        ss_res = np.sum((log_rho_m - pred)**2)
        ss_tot = np.sum((log_rho_m - log_rho_m.mean())**2)
        R2 = 1 - ss_res/ss_tot

        c0, c1, c2 = coefs
        A_pred = np.cos(theta0)
        B_pred = -(1 + np.cos(theta0))
        print(f"{label:>12s}  {c0:+8.3f} {c1:+8.4f} {c2:+8.4f}  "
              f"{A_pred:+10.4f}  {B_pred:+10.4f}  {R2:6.4f}")

    print("\nDone.")

if __name__ == "__main__":
    main()
