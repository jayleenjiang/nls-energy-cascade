#!/usr/bin/env python3
"""
Flux distribution + exponential-tail analysis
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

N_LIST = [10, 20, 30, 40]
# Fit the exponential tail over a probability window [S_LO, S_HI]:
# start below the near-Gaussian core (surv < S_HI) and stop before the far
# faster-than-exponential cutoff (surv > S_LO). This isolates the clean
# exponential middle-tail. Widen/narrow if R^2 is low or the window is empty.
S_HI = 0.20     # tail starts where survival drops below this (~80th percentile)
S_LO = 5e-3     # stop before the extreme cutoff (drop the last ~0.5%)
COLORS = ["#c0392b", "#e08a1e", "#27ae60", "#2c6fbb"]


def skew(x):
    m = x.mean(); s = x.std()
    return np.mean(((x - m) / s) ** 3) if s > 0 else 0.0


def exp_tail_fit(f, s_lo=S_LO, s_hi=S_HI):
    """Linear fit of log P[flux > x] vs x over survival in [s_lo, s_hi].
    Returns (lambda, c, R2, xs, surv, xt, fitted_survival, (xmin,xmax))."""
    xs = np.sort(f)
    surv = 1.0 - np.arange(len(xs)) / len(xs)        # P[X > xs[i]]
    m = (surv <= s_hi) & (surv >= s_lo)              # clean middle-tail window
    xt, st = xs[m], surv[m]
    A = np.vstack([xt, np.ones_like(xt)]).T
    (slope, c), *_ = np.linalg.lstsq(A, np.log(st), rcond=None)
    pred = slope * xt + c
    ss_res = np.sum((np.log(st) - pred) ** 2)
    ss_tot = np.sum((np.log(st) - np.log(st).mean()) ** 2)
    R2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return -slope, c, R2, xs, surv, xt, np.exp(pred), (xt.min(), xt.max())


def main():
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    print(f"{'n':>3} {'N':>6} {'mean':>9} {'std':>8} {'skew':>7} "
          f"{'lambda':>8} {'R2':>7}  fit-window(x)")
    rows = []
    for n, col in zip(N_LIST, COLORS):
        fn = f"flux_n{n}.txt"
        try:
            f = np.loadtxt(fn, comments="#")
        except OSError:
            print(f"  (skip {fn}: not found)")
            continue
        lam, c, R2, xs, surv, xt, fit, win = exp_tail_fit(f)
        ax[0].hist(f, bins=50, density=True, histtype="step", color=col, lw=1.6,
                   label=f"n={n} (N={len(f)})")
        ax[1].semilogy(xs, surv, "-", color=col, lw=1.2,
                       label=fr"n={n}: $\lambda$={lam:.1f}, $R^2$={R2:.4f}")
        ax[1].semilogy(xt, fit, "--", color="k", lw=0.9, alpha=0.6)
        print(f"{n:>3} {len(f):>6} {f.mean():>9.4f} {f.std():>8.4f} "
              f"{skew(f):>7.3f} {lam:>8.1f} {R2:>7.4f}  [{win[0]:.3f},{win[1]:.3f}]")
        rows.append((n, f.mean(), f.std(), lam, R2))

    ax[0].set_xlabel("time-averaged flux"); ax[0].set_ylabel("PDF")
    ax[0].set_title("Flux distribution"); ax[0].legend(fontsize=8, frameon=False)
    ax[1].set_xlabel("flux $x$"); ax[1].set_ylabel(r"$P[\,\mathrm{flux} > x\,]$")
    ax[1].set_title(fr"Exponential tail fit over survival $\in$[{S_LO:g}, {S_HI:g}]")
    ax[1].legend(fontsize=8, frameon=False)
    plt.tight_layout()
    plt.savefig("flux_distribution.png", dpi=200, bbox_inches="tight")
    print("\nWrote flux_distribution.png")

    if len(rows) >= 2:
        a = np.array(rows)
        p = np.polyfit(np.log(a[:, 0]), np.log(a[:, 1]), 1)
        print(f"\nmean flux ~ n^({p[0]:.2f})   (|exponent|~2 => diffusive/insulating)")
        pl = np.polyfit(a[:, 0], a[:, 3], 1)
        print(f"lambda ~ {pl[0]:.2f}*n + {pl[1]:.1f}")
        print(f"lambda*std = {a[:,3]*a[:,2]}  (~const => common rescaled shape)")


if __name__ == "__main__":
    main()
