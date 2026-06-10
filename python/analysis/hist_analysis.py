"""
Reader / analyzer for .hist files produced by lte_histogram.cpp.

Loads the C++-accumulated 3D histogram of (I_a, I_b, theta), reconstructs
the density, and runs the same Form A / B / C fits as lte_analysis.py.

Usage:
    python3 hist_analysis.py histo_N25_j12.hist
"""

import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt


def load_hist(fname):
    """Parse a .hist file. Returns dict with NB, ranges, counts (3D array)."""
    meta = {}
    counts = None
    with open(fname) as f:
        lines = f.readlines()

    # parse header keywords
    idx = 0
    for i, ln in enumerate(lines):
        ln = ln.strip()
        if ln.startswith("#") or ln == "":
            continue
        parts = ln.split()
        if parts[0] in ("NB",):
            meta["NB"] = int(parts[1])
        elif parts[0] in ("I_LO", "I_HI", "TH_LO", "TH_HI"):
            meta[parts[0]] = float(parts[1])
        elif parts[0] in ("TOTAL", "OVERFLOW"):
            meta[parts[0]] = int(parts[1])
        elif len(parts) == 4:
            # first data line
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
    ap.add_argument("--min-count", type=int, default=50,
                    help="drop bins below this count from the fit")
    ap.add_argument("--I-fit-max", type=float, default=3.2,
                    help="only fit bins with I below this")
    args = ap.parse_args()

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

    # bin centers
    dI = (I_HI - I_LO) / NB
    dth = (TH_HI - TH_LO) / NB
    Ia_c = I_LO + dI * (np.arange(NB) + 0.5)
    Ib_c = I_LO + dI * (np.arange(NB) + 0.5)
    th_c = TH_LO + dth * (np.arange(NB) + 0.5)

    bin_vol = dI * dI * dth
    in_count = counts.sum()
    rho = counts / (in_count * bin_vol)   # normalized over the in-range box

    IAm, IBm, THm = np.meshgrid(Ia_c, Ib_c, th_c, indexing='ij')

    # mask for fitting
    mask = (counts >= args.min_count) & \
           (IAm < args.I_fit_max) & (IBm < args.I_fit_max)
    print(f"  bins used in fit: {mask.sum():,} / {counts.size:,}")

    log_rho = np.log(rho[mask])
    a = IAm[mask]
    b = IBm[mask]
    th = THm[mask]

    ss_tot = np.sum((log_rho - log_rho.mean())**2)

    # ---- Form A ----
    F_A = -(a**2 + b**2)/4.0 + a*b*np.cos(th)
    DA = np.c_[np.ones_like(F_A), F_A]
    cA, *_ = np.linalg.lstsq(DA, log_rho, rcond=None)
    R2_A = 1 - np.sum((log_rho - DA@cA)**2)/ss_tot

    # ---- Form B ----
    F_B = (a**2 + b**2)/4.0 + a*b*(1.0 + np.cos(th))
    DB = np.c_[np.ones_like(F_B), F_B]
    cB, *_ = np.linalg.lstsq(DB, log_rho, rcond=None)
    R2_B = 1 - np.sum((log_rho - DB@cB)**2)/ss_tot

    # ---- Form C ----
    DC = np.c_[np.ones_like(a), a**2, b**2, a*b, a*b*np.cos(th)]
    cC, *_ = np.linalg.lstsq(DC, log_rho, rcond=None)
    R2_C = 1 - np.sum((log_rho - DC@cC)**2)/ss_tot

    print("\n=== FIT RESULTS ===")
    print(f"Form A: R^2 = {R2_A:.5f}   (coef = {cA[1]:+.4f})")
    print(f"Form B: R^2 = {R2_B:.5f}   (coef = {cB[1]:+.4f})")
    print(f"Form C: R^2 = {R2_C:.5f}")
    names = ['const', 'I_a^2', 'I_b^2', 'I_a*I_b', 'I_a*I_b*cos']
    for nm, c in zip(names, cC):
        print(f"   {nm:14s}: {c:+.6f}")
    print(f"   c1/c2 = {cC[1]/cC[2]:.4f}  (local Gibbs: 1.0)")
    print(f"   c3/c4 = {cC[3]/cC[4]:.4f}  (local Gibbs: 1.0)")
    print(f"   c4/c3 = {cC[4]/cC[3]:.4f}  (angle modulation strength kappa)")
    print(f"   c1/c3 = {cC[1]/cC[3]:.4f}  (local Gibbs: 0.25)")

    # ---- plot: predicted vs observed for the three forms ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    for ax, pred, R2, title in [
        (axes[0], DA@cA, R2_A, "Form A"),
        (axes[1], DB@cB, R2_B, "Form B"),
        (axes[2], DC@cC, R2_C, "Form C"),
    ]:
        ax.scatter(pred, log_rho, s=2, alpha=0.3)
        lo = min(pred.min(), log_rho.min())
        hi = max(pred.max(), log_rho.max())
        ax.plot([lo, hi], [lo, hi], 'r--', lw=1)
        ax.set_xlabel("predicted log rho")
        ax.set_ylabel("observed log rho")
        ax.set_title(f"{title}\nR^2 = {R2:.5f}")
        ax.grid(alpha=0.3)
    plt.suptitle(f"LTE fit from C++ histogram ({total:,} samples)", fontsize=12)
    plt.tight_layout()
    out = args.histfile.replace(".hist", "_fit.png")
    plt.savefig(out, dpi=120, bbox_inches='tight')
    print(f"\nSaved: {out}")

    print("\nDone.")


if __name__ == "__main__":
    main()
