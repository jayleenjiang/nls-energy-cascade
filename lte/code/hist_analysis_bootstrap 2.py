"""
Reader / analyzer for .hist files produced by lte_histogram*.cpp.

Loads the C++-accumulated 3D histogram of (I_a, I_b, theta), reconstructs
the density rho, runs Form A / B / C fits, and reports confidence
intervals using Agresti-Coull binomial CIs on the per-bin proportions.

Two reporting tracks are produced:

  TRACK 1 (pure Agresti-Coull + WLS):
      Each bin has count n_b out of N total. The Agresti-Coull
      half-width on p_b is propagated to a sigma on log rho_b, used
      as 1/sigma^2 weights in a weighted least squares fit of Form C.
      Coefficient covariance is taken from (X' W X)^{-1}.
      This is the literal "Agresti-Coull binomial CI for fit
      coefficients" interpretation.

  TRACK 2 (OLS + sandwich estimator with Agresti-Coull sigmas):
      OLS central estimates (same as Poisson bootstrap mean), but
      error bars from the Huber-White sandwich estimator
          Cov(beta_OLS) = (X'X)^{-1} (X' diag(sigma_i^2) X) (X'X)^{-1}
      with sigma_i from Agresti-Coull. Robust against Form C
      misspecification.

The two tracks generally agree on whether ratios are far from
strict-Gibbs predictions; their central values can differ by several
percent if Form C is not exactly satisfied. R^2 is reported in the
appropriate metric for each track (weighted for WLS, unweighted for
OLS). Poisson bootstrap also shown as a third cross-check.

Usage:
    python3 hist_analysis.py histo_N25_j12.hist [--n-boot 200]
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
    """Return (sigma_logp, p_hat) per bin. sigma_logp is the 1-sigma
    Agresti-Coull standard error on log p_b (half-width / z, then divided
    by p_hat via delta method d log p / dp = 1/p)."""
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("histfile")
    ap.add_argument("--min-count", type=int, default=50)
    ap.add_argument("--I-fit-max", type=float, default=3.2)
    ap.add_argument("--n-boot", type=int, default=100)
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

    # Per-bin Agresti-Coull sigma
    N_in = int(counts.sum())
    sigma_logp, p_hat = agresti_coull_sigma_logp(counts_masked, N_in, z=z95)
    log_rho_obs = np.log(p_hat) - np.log(bin_vol)

    # Design matrix for Form C
    DC = np.c_[np.ones_like(a), a**2, b**2, a*b, a*b*np.cos(th)]
    # Design matrices for Forms A, B (1-parameter, OLS shape comparison)
    F_A = -(a**2 + b**2)/4.0 + a*b*np.cos(th)
    DA = np.c_[np.ones_like(F_A), F_A]
    F_B = (a**2 + b**2)/4.0 + a*b*(1.0 + np.cos(th))
    DB = np.c_[np.ones_like(F_B), F_B]

    # ============================================================
    # Forms A and B (OLS for shape comparison only)
    # ============================================================
    cA, *_ = np.linalg.lstsq(DA, log_rho_obs, rcond=None)
    cB, *_ = np.linalg.lstsq(DB, log_rho_obs, rcond=None)
    ss_tot_unw = np.sum((log_rho_obs - log_rho_obs.mean())**2)
    R2_A_unw = 1 - np.sum((log_rho_obs - DA @ cA)**2) / ss_tot_unw
    R2_B_unw = 1 - np.sum((log_rho_obs - DB @ cB)**2) / ss_tot_unw

    # ============================================================
    # TRACK 1: Pure Agresti-Coull + WLS for Form C
    # weights = 1 / sigma_logp^2
    # ============================================================
    w = 1.0 / (sigma_logp**2)
    DTWD = (DC * w[:, None]).T @ DC
    DTWy = (DC * w[:, None]).T @ log_rho_obs
    cC_wls = np.linalg.solve(DTWD, DTWy)
    cov_wls = np.linalg.inv(DTWD)
    sig_wls = np.sqrt(np.diag(cov_wls))

    # WLS-appropriate R^2: weighted
    wmean = np.sum(w * log_rho_obs) / np.sum(w)
    ss_tot_w = np.sum(w * (log_rho_obs - wmean)**2)
    ss_res_w = np.sum(w * (log_rho_obs - DC @ cC_wls)**2)
    R2_wls_weighted = 1 - ss_res_w / ss_tot_w

    # WLS coefficients evaluated under unweighted R^2 (for comparison)
    R2_wls_unw = 1 - np.sum((log_rho_obs - DC @ cC_wls)**2) / ss_tot_unw

    # ============================================================
    # TRACK 2: OLS + Agresti-Coull sandwich
    # ============================================================
    DtD = DC.T @ DC
    DtD_inv = np.linalg.inv(DtD)
    cC_ols = DtD_inv @ (DC.T @ log_rho_obs)
    Omega_diag = sigma_logp**2
    middle = (DC.T * Omega_diag) @ DC
    cov_ols_sw = DtD_inv @ middle @ DtD_inv
    sig_ols_sw = np.sqrt(np.diag(cov_ols_sw))

    R2_ols_unw = 1 - np.sum((log_rho_obs - DC @ cC_ols)**2) / ss_tot_unw
    R2_ols_w   = 1 - np.sum(w * (log_rho_obs - DC @ cC_ols)**2) / ss_tot_w

    # ============================================================
    # Print
    # ============================================================
    print("\n=== Forms A and B (OLS, shape comparison) ===")
    print(f"Form A: R^2 (unweighted) = {R2_A_unw:.5f}")
    print(f"Form B: R^2 (unweighted) = {R2_B_unw:.5f}")

    print("\n=== TRACK 1: Pure Agresti-Coull + WLS (Form C) ===")
    print(f"  R^2 weighted   = {R2_wls_weighted:.5f}   (the metric WLS optimizes)")
    print(f"  R^2 unweighted = {R2_wls_unw:.5f}   (not the WLS criterion)")
    print(f"  c0 = {cC_wls[0]:+.5f} +- {sig_wls[0]:.5f}")
    print(f"  c1 = {cC_wls[1]:+.5f} +- {sig_wls[1]:.5f}")
    print(f"  c2 = {cC_wls[2]:+.5f} +- {sig_wls[2]:.5f}")
    print(f"  c3 = {cC_wls[3]:+.5f} +- {sig_wls[3]:.5f}")
    print(f"  c4 = {cC_wls[4]:+.5f} +- {sig_wls[4]:.5f}")
    r12, s12 = delta_ratio(cC_wls[1], cC_wls[2],
                           cov_wls[1,1], cov_wls[2,2], cov_wls[1,2])
    r34, s34 = delta_ratio(cC_wls[3], cC_wls[4],
                           cov_wls[3,3], cov_wls[4,4], cov_wls[3,4])
    r13, s13 = delta_ratio(cC_wls[1], cC_wls[3],
                           cov_wls[1,1], cov_wls[3,3], cov_wls[1,3])
    print(f"  c1/c2 = {r12:.5f} +- {s12:.5f}   95% CI=[{r12-z95*s12:.5f}, {r12+z95*s12:.5f}]")
    print(f"  c3/c4 = {r34:.5f} +- {s34:.5f}   95% CI=[{r34-z95*s34:.5f}, {r34+z95*s34:.5f}]")
    print(f"  c1/c3 = {r13:.5f} +- {s13:.5f}   95% CI=[{r13-z95*s13:.5f}, {r13+z95*s13:.5f}]")
    print(f"  Strict Gibbs targets: c1/c2=1, c3/c4=1, c1/c3=0.25")
    print(f"  z vs strict: c1/c2={(r12-1)/s12:+.1f}  c3/c4={(r34-1)/s34:+.1f}  "
          f"c1/c3={(r13-0.25)/s13:+.1f}")

    print("\n=== TRACK 2: OLS + Agresti-Coull sandwich (Form C) ===")
    print(f"  R^2 unweighted = {R2_ols_unw:.5f}   (the metric OLS optimizes)")
    print(f"  R^2 weighted   = {R2_ols_w:.5f}")
    print(f"  c0 = {cC_ols[0]:+.5f} +- {sig_ols_sw[0]:.5f}")
    print(f"  c1 = {cC_ols[1]:+.5f} +- {sig_ols_sw[1]:.5f}")
    print(f"  c2 = {cC_ols[2]:+.5f} +- {sig_ols_sw[2]:.5f}")
    print(f"  c3 = {cC_ols[3]:+.5f} +- {sig_ols_sw[3]:.5f}")
    print(f"  c4 = {cC_ols[4]:+.5f} +- {sig_ols_sw[4]:.5f}")
    r12o, s12o = delta_ratio(cC_ols[1], cC_ols[2],
                             cov_ols_sw[1,1], cov_ols_sw[2,2], cov_ols_sw[1,2])
    r34o, s34o = delta_ratio(cC_ols[3], cC_ols[4],
                             cov_ols_sw[3,3], cov_ols_sw[4,4], cov_ols_sw[3,4])
    r13o, s13o = delta_ratio(cC_ols[1], cC_ols[3],
                             cov_ols_sw[1,1], cov_ols_sw[3,3], cov_ols_sw[1,3])
    print(f"  c1/c2 = {r12o:.5f} +- {s12o:.5f}   95% CI=[{r12o-z95*s12o:.5f}, {r12o+z95*s12o:.5f}]")
    print(f"  c3/c4 = {r34o:.5f} +- {s34o:.5f}   95% CI=[{r34o-z95*s34o:.5f}, {r34o+z95*s34o:.5f}]")
    print(f"  c1/c3 = {r13o:.5f} +- {s13o:.5f}   95% CI=[{r13o-z95*s13o:.5f}, {r13o+z95*s13o:.5f}]")
    print(f"  z vs strict: c1/c2={(r12o-1)/s12o:+.1f}  c3/c4={(r34o-1)/s34o:+.1f}  "
          f"c1/c3={(r13o-0.25)/s13o:+.1f}")

    # ============================================================
    # Track 3 (cross-check): Poisson bootstrap, OLS-style
    # ============================================================
    print(f"\n=== Cross-check: Poisson bootstrap ({args.n_boot} resamples, OLS each) ===")
    rng = np.random.default_rng(12345)
    boot_c34 = np.empty(args.n_boot)
    boot_c13 = np.empty(args.n_boot)
    boot_c12 = np.empty(args.n_boot)
    for k in range(args.n_boot):
        c_resamp = rng.poisson(counts_masked).astype(float)
        keep = c_resamp > 0
        if keep.sum() < 100:
            boot_c34[k] = boot_c13[k] = boot_c12[k] = np.nan; continue
        ck = c_resamp[keep]
        ak = a[keep]; bk = b[keep]; thk = th[keep]
        rho_k = ck / (ck.sum() * bin_vol)
        log_rho_k = np.log(rho_k)
        Dk = np.c_[np.ones_like(ak), ak**2, bk**2, ak*bk, ak*bk*np.cos(thk)]
        ck_C, *_ = np.linalg.lstsq(Dk, log_rho_k, rcond=None)
        boot_c12[k] = ck_C[1]/ck_C[2]
        boot_c34[k] = ck_C[3]/ck_C[4]
        boot_c13[k] = ck_C[1]/ck_C[3]
    def br(label, arr):
        m=np.nanmean(arr); s=np.nanstd(arr,ddof=1)
        lo,hi=np.nanpercentile(arr,[2.5,97.5])
        print(f"  {label}: boot mean={m:+.5f}  std={s:.5f}  95% pct=[{lo:.5f}, {hi:.5f}]")
    br("c1/c2", boot_c12)
    br("c3/c4", boot_c34)
    br("c1/c3", boot_c13)

    print("\nNotes:")
    print("  WLS uses weights = 1/sigma_AC^2 and is the literal Agresti-Coull-")
    print("  weighted fit; its R^2 should be read on the weighted metric.  Its")
    print("  central estimates can differ from OLS / bootstrap because the AC")
    print("  weights ~ bin count, biasing toward the high-density corner.")
    print("  OLS + sandwich gives OLS central values (matching bootstrap mean)")
    print("  with heteroscedastic-consistent error bars built from the AC sigmas.")
    print("  All methods ignore trajectory autocorrelation.")

    # ============================================================
    # Plot: side-by-side WLS and OLS ratios
    # ============================================================
    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    def panel(ax, r, s, target, label, color, title):
        x = np.linspace(r - 4*s - 1e-6, r + 4*s + 1e-6, 400)
        y = np.exp(-0.5*((x-r)/s)**2) / (s*np.sqrt(2*np.pi)) if s>0 else np.zeros_like(x)
        ax.plot(x, y, color=color, lw=2)
        if s > 0:
            ax.fill_between(x, 0, y, where=(x>=r-z95*s)&(x<=r+z95*s),
                            color=color, alpha=0.3, label='95% CI')
        ax.axvline(target, color='r', ls='--', label='strict Gibbs')
        ax.axvline(r, color='k', lw=1.5, label='central')
        ax.set_xlabel(label); ax.legend(fontsize=8); ax.grid(alpha=0.3)
        ax.set_title(title, fontsize=10)

    panel(axes[0,0], r12, s12, 1.0,  "c1/c2", 'C0', "WLS  c1/c2")
    panel(axes[0,1], r34, s34, 1.0,  "c3/c4", 'C1', "WLS  c3/c4")
    panel(axes[0,2], r13, s13, 0.25, "c1/c3", 'C2', "WLS  c1/c3")
    panel(axes[1,0], r12o, s12o, 1.0,  "c1/c2", 'C0', "OLS+sandwich  c1/c2")
    panel(axes[1,1], r34o, s34o, 1.0,  "c3/c4", 'C1', "OLS+sandwich  c3/c4")
    panel(axes[1,2], r13o, s13o, 0.25, "c1/c3", 'C2', "OLS+sandwich  c1/c3")

    plt.suptitle(f"AC-WLS (top) vs OLS+AC-sandwich (bottom) -- {args.histfile} ({total:,} samples)",
                 fontsize=12)
    plt.tight_layout()
    out = args.histfile.replace(".hist", "_fit.png")
    plt.savefig(out, dpi=120, bbox_inches='tight')
    print(f"\nSaved: {out}")
    print("Done.")


if __name__ == "__main__":
    main()
