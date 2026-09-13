"""
Combined analysis on full-chain trajectory CSV.

Reads traj_N25_full.csv (must have all modes 0..n-1 dumped) and produces:

  1. Thermal profile <I_j> vs j across the chain
  2. Multi-site LTE scan: fit Form C at each interior site j,
     report the angle-coupling ratio c_3/c_4 vs j.

Usage:
    python lte_full_analysis.py traj_N25_full.csv [--n 25]
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def fit_form_C(Ia, Ib, theta, bins=40, min_count=20, Imax=None):
    """Fit log rho = const + c1*Ia^2 + c2*Ib^2 + c3*Ia*Ib + c4*Ia*Ib*cos(theta)
    Returns (coefs[5], R^2, n_bins_used)."""
    if Imax is None:
        Imax = np.percentile(np.r_[Ia, Ib], 99)
    Imin = 1e-3
    keep = (Ia > Imin) & (Ib > Imin) & (Ia < Imax) & (Ib < Imax)
    Ia, Ib, theta = Ia[keep], Ib[keep], theta[keep]

    Ia_edges = np.linspace(Imin, Imax, bins+1)
    Ib_edges = np.linspace(Imin, Imax, bins+1)
    th_edges = np.linspace(-np.pi, np.pi, bins+1)
    H, _ = np.histogramdd(np.c_[Ia, Ib, theta],
                          bins=(Ia_edges, Ib_edges, th_edges))
    Ia_c = 0.5*(Ia_edges[1:] + Ia_edges[:-1])
    Ib_c = 0.5*(Ib_edges[1:] + Ib_edges[:-1])
    th_c = 0.5*(th_edges[1:] + th_edges[:-1])
    dvol = (Ia_edges[1]-Ia_edges[0]) * (Ib_edges[1]-Ib_edges[0]) * (th_edges[1]-th_edges[0])
    rho = H / (H.sum() * dvol)
    IAm, IBm, THm = np.meshgrid(Ia_c, Ib_c, th_c, indexing='ij')

    mask = H >= min_count
    log_rho = np.log(rho[mask])
    Iam, Ibm, thm = IAm[mask], IBm[mask], THm[mask]

    feats = np.c_[
        np.ones_like(Iam),
        Iam**2,
        Ibm**2,
        Iam*Ibm,
        Iam*Ibm*np.cos(thm),
    ]
    coefs, *_ = np.linalg.lstsq(feats, log_rho, rcond=None)
    pred = feats @ coefs
    ss_res = np.sum((log_rho - pred)**2)
    ss_tot = np.sum((log_rho - log_rho.mean())**2)
    R2 = 1 - ss_res/ss_tot
    return coefs, R2, mask.sum()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--n", type=int, default=25, help="chain length")
    ap.add_argument("--bins", type=int, default=40)
    ap.add_argument("--min-count", type=int, default=20)
    args = ap.parse_args()

    n = args.n
    print(f"Loading {args.csv} ...")
    df = pd.read_csv(args.csv)
    print(f"  {len(df)} samples")

    # check columns
    needed = [f"I_{j}" for j in range(n)] + [f"phi_{j}" for j in range(n)]
    missing = [c for c in needed if c not in df.columns]
    if missing:
        print(f"ERROR: missing columns {missing[:5]}...")
        print("Did you run with j_lo=0 j_hi={}?".format(n-1))
        return

    # ============================================================
    # 1. Thermal profile <I_j>
    # ============================================================
    print("\n=== THERMAL PROFILE ===")
    mean_I = np.array([df[f"I_{j}"].mean() for j in range(n)])
    std_I  = np.array([df[f"I_{j}"].std()  for j in range(n)])

    print(f"{'j':>4s} {'<I_j>':>12s} {'std':>12s}")
    for j in range(n):
        print(f"{j:4d} {mean_I[j]:12.4f} {std_I[j]:12.4f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].plot(range(n), mean_I, 'o-', lw=1.5)
    axes[0].set_xlabel("mode index j")
    axes[0].set_ylabel("<I_j>")
    axes[0].set_title(f"Thermal profile (T_1=10, T_n=2, n={n})")
    axes[0].grid(alpha=0.3)

    axes[1].semilogy(range(n), mean_I, 'o-', lw=1.5)
    axes[1].set_xlabel("mode index j")
    axes[1].set_ylabel("<I_j> (log scale)")
    axes[1].set_title("Thermal profile (log y)")
    axes[1].grid(alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig("thermal_profile.png", dpi=120, bbox_inches='tight')
    print("Saved: thermal_profile.png")

    # ============================================================
    # 2. Multi-site LTE scan (Form C fit)
    # ============================================================
    print("\n=== MULTI-SITE LTE SCAN (Form C) ===")
    print("Fitting log rho ~ const + c1 I_a^2 + c2 I_b^2 + c3 I_a I_b + c4 I_a I_b cos(theta)")
    print(f"\n{'j':>3s} {'R^2':>8s} {'c1':>10s} {'c2':>10s} {'c3':>10s} "
          f"{'c4':>10s} {'c1/c2':>8s} {'c3/c4':>8s} {'c1/c3':>8s} {'<I_j>':>8s}")

    results = []
    for j in range(n-1):  # need pair (j, j+1)
        Ia = df[f"I_{j}"].to_numpy()
        Ib = df[f"I_{j+1}"].to_numpy()
        pa = df[f"phi_{j}"].to_numpy()
        pb = df[f"phi_{j+1}"].to_numpy()
        theta = 2.0*(pb - pa)
        theta = np.mod(theta + np.pi, 2*np.pi) - np.pi
        try:
            coefs, R2, nbins = fit_form_C(Ia, Ib, theta,
                                          bins=args.bins,
                                          min_count=args.min_count)
            c0, c1, c2, c3, c4 = coefs
            r12 = c1/c2 if abs(c2) > 1e-9 else np.nan
            r34 = c3/c4 if abs(c4) > 1e-9 else np.nan
            r13 = c1/c3 if abs(c3) > 1e-9 else np.nan
            print(f"{j:3d} {R2:8.4f} {c1:+10.4f} {c2:+10.4f} {c3:+10.4f} "
                  f"{c4:+10.4f} {r12:8.3f} {r34:8.3f} {r13:8.3f} "
                  f"{mean_I[j]:8.3f}")
            results.append((j, R2, c0, c1, c2, c3, c4))
        except Exception as e:
            print(f"{j:3d} FAILED: {e}")
            results.append((j, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan))

    res = np.array(results)
    js = res[:, 0]
    R2s = res[:, 1]
    c1s, c2s, c3s, c4s = res[:, 3], res[:, 4], res[:, 5], res[:, 6]

    # ============================================================
    # 3. Plot LTE diagnostics
    # ============================================================
    fig, axes = plt.subplots(2, 2, figsize=(14, 9))

    # (a) R^2 vs site
    axes[0, 0].plot(js, R2s, 'o-', lw=1.5)
    axes[0, 0].set_xlabel("pair index j (modes j, j+1)")
    axes[0, 0].set_ylabel("R^2 (Form C fit)")
    axes[0, 0].set_title("LTE fit quality across chain")
    axes[0, 0].grid(alpha=0.3)
    axes[0, 0].set_ylim(0, 1.05)

    # (b) c1/c2 (symmetry between I_a, I_b)
    axes[0, 1].plot(js, c1s/c2s, 'o-', lw=1.5)
    axes[0, 1].axhline(1.0, color='r', ls='--', label='strict symmetry (c1=c2)')
    axes[0, 1].set_xlabel("pair index j")
    axes[0, 1].set_ylabel("c1 / c2  (= I_a^2 / I_b^2 coef ratio)")
    axes[0, 1].set_title("Symmetry between I_j, I_{j+1}")
    axes[0, 1].legend()
    axes[0, 1].grid(alpha=0.3)

    # (c) c3/c4 (angle coupling strength relative to no-cos cross term)
    axes[1, 0].plot(js, c3s/c4s, 'o-', lw=1.5, label='c3/c4 (data)')
    axes[1, 0].axhline(1.0, color='r', ls='--', label='strict local Gibbs (c3=c4)')
    axes[1, 0].set_xlabel("pair index j")
    axes[1, 0].set_ylabel("c3 / c4")
    axes[1, 0].set_title("Angle coupling: c3/c4=1 means strict local Gibbs;\n"
                         ">1 means cos(theta) modulation suppressed")
    axes[1, 0].legend()
    axes[1, 0].grid(alpha=0.3)

    # (d) implied 'temperature' from c3 (= -beta in local Gibbs interpretation)
    Tloc = -1.0/(2.0*c3s) * 4.0  # if c3 = -beta, T = 1/beta; but our parametrization has beta in c1=beta/4 and c3=beta/1 => use whichever
    # actually for Form B: log rho = const + b * F_B,  b = -beta
    # in Form C, c1 = -beta/4 and c3 = -beta. So beta = -c3 = -4*c1.
    # T = 1/beta. Use c3 estimate.
    beta_from_c3 = -c3s
    T_from_c3 = 1.0/beta_from_c3
    # linear interpolation profile (T_1=10 at j=0, T_n=2 at j=n-1)
    j_axis = np.arange(n)
    T_lin = 10.0 - (10.0-2.0) * j_axis / (n-1)
    axes[1, 1].plot(js, T_from_c3, 'o-', lw=1.5, label='T_loc from c3 = -beta')
    axes[1, 1].plot(j_axis, T_lin, 'r--', label='linear T_1->T_n')
    axes[1, 1].set_xlabel("pair index j")
    axes[1, 1].set_ylabel("implied local temperature")
    axes[1, 1].set_title("Local 'temperature' from fit (if local Gibbs holds)")
    axes[1, 1].legend()
    axes[1, 1].grid(alpha=0.3)

    plt.tight_layout()
    plt.savefig("lte_scan.png", dpi=120, bbox_inches='tight')
    print("\nSaved: lte_scan.png")

    # save numeric results
    out = pd.DataFrame({
        'j': js.astype(int),
        'R2': R2s,
        'c0_const': res[:, 2],
        'c1_Ia2':   c1s,
        'c2_Ib2':   c2s,
        'c3_IaIb':  c3s,
        'c4_IaIbcos': c4s,
        'c1_over_c2': c1s/c2s,
        'c3_over_c4': c3s/c4s,
        'c1_over_c3': c1s/c3s,
        'mean_I_j':  mean_I[:n-1],
    })
    out.to_csv("lte_scan_results.csv", index=False)
    print("Saved: lte_scan_results.csv")

    print("\nDone.")

if __name__ == "__main__":
    main()
