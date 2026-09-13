"""
Agresti-Coull binomial confidence intervals on per-bin proportions, propagated to a weighted least-squares fit of Form C.

For each bin b with count n_b out of N total samples in the fit region,
treat n_b as Binomial(N, p_b). The Agresti-Coull 95% CI on p_b is
    p~_b      = (n_b + 2) / (N + 4)
    halfwidth = z * sqrt(p~_b (1 - p~_b) / (N + 4)),   z = 1.96
The 1-sigma standard error on log p_b is propagated by delta method:
    sigma(log p_b) = (halfwidth / z) / p_hat_b
These sigmas are used as 1/sigma^2 weights in a WLS fit of
    log rho(I_a, I_b, theta)
       = c0 + c1 I_a^2 + c2 I_b^2 + c3 I_a I_b + c4 I_a I_b cos(theta).
Coefficient covariance is taken from (X^T W X)^{-1}; ratios c_i/c_j and
their 1-sigma errors come from the delta method.

Forms A and B (1-parameter shape checks) are also fit by WLS with the
same Agresti-Coull weights, for R^2 comparison on the same metric.

Usage:
    python3 hist_analysis.py histo_N25_j12.hist
"""

import argparse
import numpy as np
import matplotlib.pyplot as plt


def load_hist(fname):
    """Parse a .hist file. Returns dict with NB, ranges, counts (3D array)."""
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


def agresti_coull_sigma_logp(n_bin, N_total, z=1.96):
    """Agresti-Coull 1-sigma standard error on log p_b for each bin.

    Returns (sigma_logp, p_hat) where
        sigma_logp[b] = (z * sqrt(p~(1-p~)/n~) / z) / p_hat
                      = sqrt(p~(1-p~)/n~) / p_hat
    is the 1-sigma equivalent (95% half-width divided by z, then divided
    by p_hat to convert from sigma on p to sigma on log p).
    """
    n_tilde = N_total + 4.0
    p_tilde = (n_bin + 2.0) / n_tilde
    halfwidth_p = z * np.sqrt(p_tilde * (1.0 - p_tilde) / n_tilde)
    p_hat = n_bin.astype(float) / N_total
    safe_p = np.where(p_hat > 0, p_hat, 1.0)
    sigma_logp = (halfwidth_p / z) / safe_p
    return sigma_logp, p_hat


def delta_ratio(x, y, Vxx, Vyy, Vxy):
    """Var(x/y) by the delta method."""
    var = Vxx / y**2 + (x**2 * Vyy) / y**4 - (2.0 * x * Vxy) / y**3
    return x / y, np.sqrt(max(var, 0.0))


def wls_fit(D, y, w):
    """Weighted least squares with weights w (= 1/sigma^2).
    Returns (beta, cov_beta)."""
    DTWD = (D * w[:, None]).T @ D
    DTWy = (D * w[:, None]).T @ y
    beta = np.linalg.solve(DTWD, DTWy)
    cov = np.linalg.inv(DTWD)
    return beta, cov


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("histfile")
    ap.add_argument("--min-count", type=int, default=50,
                    help="drop bins below this count from the fit")
    ap.add_argument("--I-fit-max", type=float, default=3.2,
                    help="only fit bins with I below this")
    args = ap.parse_args()
    z95 = 1.96

    print(f"Loading {args.histfile} ...")
    H = load_hist(args.histfile)
    NB = H["NB"]
    counts = H["counts"]
    I_LO, I_HI = H["I_LO"], H["I_HI"]
    TH_LO, TH_HI = H["TH_LO"], H["TH_HI"]
    total = H.get("TOTAL", counts.sum())
    overflow = H.get("OVERFLOW", 0)
    print(f"  NB={NB}, total={total:,}, overflow={overflow:,} "
          f"({100*overflow/max(total,1):.2f}%)")

    dI = (I_HI - I_LO) / NB
    dth = (TH_HI - TH_LO) / NB
    Ia_c = I_LO + dI * (np.arange(NB) + 0.5)
    Ib_c = I_LO + dI * (np.arange(NB) + 0.5)
    th_c = TH_LO + dth * (np.arange(NB) + 0.5)
    bin_vol = dI * dI * dth

    IAm, IBm, THm = np.meshgrid(Ia_c, Ib_c, th_c, indexing='ij')
    mask = (counts >= args.min_count) & \
           (IAm < args.I_fit_max) & (IBm < args.I_fit_max)
    print(f"  bins used in fit: {mask.sum():,} / {counts.size:,}")

    counts_masked = counts[mask]
    a = IAm[mask].astype(float)
    b = IBm[mask].astype(float)
    th = THm[mask].astype(float)

    # Per-bin Agresti-Coull sigma on log p
    N_in = int(counts.sum())
    sigma_logp, p_hat = agresti_coull_sigma_logp(counts_masked, N_in, z=z95)
    log_rho_obs = np.log(p_hat) - np.log(bin_vol)
    w = 1.0 / sigma_logp**2

    # Weighted total-sum-of-squares (the metric WLS optimizes against)
    wmean = np.sum(w * log_rho_obs) / np.sum(w)
    ss_tot_w = np.sum(w * (log_rho_obs - wmean)**2)

    # ============================================================
    # Form A (1-parameter): conjecture
    # log rho = a0 + a1 * [ -(I_a^2+I_b^2)/4 + I_a I_b cos(theta) ]
    # ============================================================
    F_A = -(a**2 + b**2)/4.0 + a*b*np.cos(th)
    DA = np.c_[np.ones_like(F_A), F_A]
    cA, covA = wls_fit(DA, log_rho_obs, w)
    R2_A_w = 1 - np.sum(w * (log_rho_obs - DA @ cA)**2) / ss_tot_w

    # ============================================================
    # Form B (1-parameter): strict local Gibbs
    # log rho = b0 + b1 * [ (I_a^2+I_b^2)/4 + I_a I_b (1 + cos(theta)) ]
    # ============================================================
    F_B = (a**2 + b**2)/4.0 + a*b*(1.0 + np.cos(th))
    DB = np.c_[np.ones_like(F_B), F_B]
    cB, covB = wls_fit(DB, log_rho_obs, w)
    R2_B_w = 1 - np.sum(w * (log_rho_obs - DB @ cB)**2) / ss_tot_w

    # ============================================================
    # Form C (4-parameter): general quadratic
    # log rho = c0 + c1 I_a^2 + c2 I_b^2 + c3 I_a I_b + c4 I_a I_b cos(theta)
    # ============================================================
    DC = np.c_[np.ones_like(a), a**2, b**2, a*b, a*b*np.cos(th)]
    cC, covC = wls_fit(DC, log_rho_obs, w)
    sigC = np.sqrt(np.diag(covC))
    R2_C_w = 1 - np.sum(w * (log_rho_obs - DC @ cC)**2) / ss_tot_w

    # Ratios + delta-method standard errors
    r12, s12 = delta_ratio(cC[1], cC[2], covC[1,1], covC[2,2], covC[1,2])
    r34, s34 = delta_ratio(cC[3], cC[4], covC[3,3], covC[4,4], covC[3,4])
    r13, s13 = delta_ratio(cC[1], cC[3], covC[1,1], covC[3,3], covC[1,3])

    # ============================================================
    # Print
    # ============================================================
    print("\n=== Form A (conjecture), WLS with Agresti-Coull weights ===")
    print(f"  R^2 (weighted) = {R2_A_w:.5f}")

    print("\n=== Form B (strict local Gibbs), WLS with Agresti-Coull weights ===")
    print(f"  R^2 (weighted) = {R2_B_w:.5f}")

    print("\n=== Form C (general quadratic), WLS with Agresti-Coull weights ===")
    print(f"  R^2 (weighted) = {R2_C_w:.5f}")
    print(f"  c0 (const)       = {cC[0]:+.5f} +- {sigC[0]:.5f}")
    print(f"  c1 (I_a^2)       = {cC[1]:+.5f} +- {sigC[1]:.5f}")
    print(f"  c2 (I_b^2)       = {cC[2]:+.5f} +- {sigC[2]:.5f}")
    print(f"  c3 (I_a I_b)     = {cC[3]:+.5f} +- {sigC[3]:.5f}")
    print(f"  c4 (I_a I_b cos) = {cC[4]:+.5f} +- {sigC[4]:.5f}")

    print("\n=== Ratios with 95% CI (delta method on WLS covariance) ===")
    def ci(r, s): return r - z95*s, r + z95*s
    lo, hi = ci(r12, s12)
    print(f"  c1/c2 = {r12:.5f} +- {s12:.5f}   95% CI = [{lo:.5f}, {hi:.5f}]   "
          f"strict Gibbs: 1.00")
    lo, hi = ci(r34, s34)
    print(f"  c3/c4 = {r34:.5f} +- {s34:.5f}   95% CI = [{lo:.5f}, {hi:.5f}]   "
          f"strict Gibbs: 1.00")
    lo, hi = ci(r13, s13)
    print(f"  c1/c3 = {r13:.5f} +- {s13:.5f}   95% CI = [{lo:.5f}, {hi:.5f}]   "
          f"strict Gibbs: 0.25")

    print("\nDistance from strict local Gibbs:")
    print(f"  c1/c2:  z = {(r12 - 1.0)/s12:+.1f}    ({100*(r12-1.0):+.2f}% off)")
    print(f"  c3/c4:  z = {(r34 - 1.0)/s34:+.1f}    ({100*(r34-1.0):+.2f}% off)")
    print(f"  c1/c3:  z = {(r13 - 0.25)/s13:+.1f}    ({100*(r13-0.25)/0.25:+.2f}% off)")

    # ============================================================
    # Plot
    # ============================================================
    fig, axes = plt.subplots(1, 4, figsize=(20, 4.5))

    # Predicted vs observed for Form C
    pred_C = DC @ cC
    axes[0].scatter(pred_C, log_rho_obs, s=2, alpha=0.3)
    lo_y, hi_y = log_rho_obs.min(), log_rho_obs.max()
    axes[0].plot([lo_y, hi_y], [lo_y, hi_y], 'r--', lw=1)
    axes[0].set_xlabel("predicted log rho (Form C, WLS)")
    axes[0].set_ylabel("observed log rho")
    axes[0].set_title(f"Form C fit  R^2_w = {R2_C_w:.4f}")
    axes[0].grid(alpha=0.3)

    def panel_ratio(ax, r, s, target, label, color):
        x = np.linspace(r - 4*s - 1e-9, r + 4*s + 1e-9, 400)
        y = np.exp(-0.5*((x-r)/s)**2) / (s*np.sqrt(2*np.pi)) if s > 0 else np.zeros_like(x)
        ax.plot(x, y, color=color, lw=2)
        if s > 0:
            ax.fill_between(x, 0, y, where=(x >= r-z95*s) & (x <= r+z95*s),
                            color=color, alpha=0.3, label='95% CI')
        ax.axvline(target, color='r', ls='--', label='strict local Gibbs')
        ax.axvline(r, color='k', lw=1.5, label='central')
        ax.set_xlabel(label); ax.legend(fontsize=8); ax.grid(alpha=0.3)

    panel_ratio(axes[1], r12, s12, 1.0,  "c1 / c2", 'C0')
    axes[1].set_title("c1/c2 (AC-WLS)")
    panel_ratio(axes[2], r34, s34, 1.0,  "c3 / c4", 'C1')
    axes[2].set_title("c3/c4 (AC-WLS)")
    panel_ratio(axes[3], r13, s13, 0.25, "c1 / c3", 'C2')
    axes[3].set_title("c1/c3 (AC-WLS)")

    plt.suptitle(f"Agresti-Coull + WLS  --  {args.histfile} ({total:,} samples)",
                 fontsize=12)
    plt.tight_layout()
    out = args.histfile.replace(".hist", "_fit.png")
    plt.savefig(out, dpi=120, bbox_inches='tight')
    print(f"\nSaved: {out}")
    print("Done.")


if __name__ == "__main__":
    main()
