"""4.3 extension, A3: linear theory of the n = 4 middle-pair trap (paper normalisation, canonical bath BC1).

Model (same as nls_cartesian_v2 / nls_boundary_growth_v2):  dc_j = (i - gamma*[j end]) F_j dt + sqrt(2 gamma T) dW,
F_j = dH/d conj(c_j) = (2M - I_j) c_j + 2 (c_{j-1}^2 + c_{j+1}^2) conj(c_j).

Pair = (c2, c3) with endpoints ~ 0.  Endpoint 1 obeys, to linear order,
    dc1 = (i - gamma) (2M c1 + 2 c2^2 conj(c1)) dt + sqrt(2 gamma T1) dW.
(1) Frozen-pair criterion.  With c2 = sqrt(I2) e^{i phi2}, phi2' = w2 = 2M - I2 + 2 I3 cos th3, u = c1 e^{-i phi2}:
    u' = i Delta u + i kappa conj(u)  (+ damping),  Delta = 2M - w2 = I2 - 2 I3 cos th3,  kappa = 2 I2.
    Endpoint 1 grows iff |Delta| < kappa, with d log I1/dt = 2 sqrt(kappa^2 - Delta^2); else it oscillates.
    Single hot mode (I3 = 0): Delta = I2, kappa = 2 I2: rate 2 sqrt(3) I2, angle 2pi/3.
    Equal anti-phase pair: Delta = 3I, kappa = 2I: stable.
(2) Exact version: Lyapunov exponent of the endpoint along the deterministic dimer orbit (gamma = 0), on a grid of
    initial (r = I3/I2, th3), I2 + I3 = 2.
(3) Leak of a stable trapped state: dM/dt = sum_ends 2 gamma (2T - I dH/dI) exactly (Ito; H conserves M).  For a
    frozen stable pair and gamma -> 0 the endpoint law is Gaussian with covariance fixed by the quadratic form
    Q = (Delta + kappa) x^2 + (Delta - kappa) y^2, giving per end
        F_end = 2 gamma T f,   f = kappa^2 (2M - Delta) / (M (Delta^2 - kappa^2)),
    i.e. f = 0.4 for the equal anti-phase pair: dM/dt = -0.8 gamma (T1 + Tn), independent of M.
    Checked here against the exact Lyapunov equation at finite gamma.
Usage: python3 theory_trap.py OUTDIR
"""
import json, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.linalg import solve_continuous_lyapunov
O = sys.argv[1]
res = {}


def frozen(I2, I3, c3):
    M = I2 + I3
    D = I2 - 2 * I3 * c3
    k = 2 * I2
    grow = 2 * np.sqrt(np.clip(k ** 2 - D ** 2, 0, None))
    with np.errstate(divide="ignore", invalid="ignore"):
        f = np.where(np.abs(D) > k, k ** 2 * (2 * M - D) / (M * (D ** 2 - k ** 2)), np.nan)
    return D, k, grow, f


def f_exact(I2, I3, c3, gamma):
    """stationary endpoint covariance at finite gamma (frozen pair), f = flux/(2 gamma T)"""
    M = I2 + I3; D = I2 - 2 * I3 * c3; k = 2 * I2
    A = np.array([[-gamma * (2 * M + k), -(D - k)], [D + k, -gamma * (2 * M - k)]])
    S = solve_continuous_lyapunov(A, -2 * gamma * np.eye(2))     # T = 1
    if np.max(np.linalg.eigvals(A).real) >= 0:
        return np.nan
    E = (2 * M + k) * S[0, 0] + (2 * M - k) * S[1, 1]
    return float(E - 2), float(S[0, 0] + S[1, 1])


# (3) leak constant, gamma dependence
res["f_equal_pair_gamma0"] = float(frozen(1.0, 1.0, -1.0)[3])
res["f_equal_pair_exact"] = {str(g): f_exact(1.0, 1.0, -1.0, g)[0] for g in [0.01, 0.05, 0.1, 0.2, 0.5]}
res["endpoint_mean_action_equal_pair_times_I_over_T"] = f_exact(1.0, 1.0, -1.0, 0.1)[1]   # <I1> = this * T / I

# (2) exact endpoint Lyapunov exponent along dimer orbits
GAMMAS = np.array([0.0, 0.02, 0.05, 0.1])[:, None]       # endpoint damping (bath) for the tangent vectors


def rhs(c2, c3, e1, e4):
    """e1, e4: shape (len(GAMMAS), grid): linearised endpoints with damping (i - gamma)"""
    I2, I3 = abs(c2) ** 2, abs(c3) ** 2
    M = I2 + I3
    d2 = 1j * ((2 * M - I2) * c2 + 2 * c3 ** 2 * np.conj(c2))
    d3 = 1j * ((2 * M - I3) * c3 + 2 * c2 ** 2 * np.conj(c3))
    d1 = (1j - GAMMAS) * (2 * M * e1 + 2 * c2 ** 2 * np.conj(e1))
    d4 = (1j - GAMMAS) * (2 * M * e4 + 2 * c3 ** 2 * np.conj(e4))
    return d2, d3, d1, d4


nr, nt = 61, 72
rr = np.exp(np.linspace(np.log(1 / 6), np.log(6), nr))
tt = np.linspace(0, 2 * np.pi, nt, endpoint=False) + np.pi / nt
R, TH = np.meshgrid(rr, tt, indexing="ij")
I2 = 2 / (1 + R); I3 = 2 - I2
c2 = np.sqrt(I2).astype(complex).ravel()
c3 = (np.sqrt(I3) * np.exp(0.5j * TH)).ravel()            # th3 = 2(phi3 - phi2)
rng = np.random.default_rng(5)
e1 = np.exp(2j * np.pi * rng.random((len(GAMMAS), c2.size))); e4 = np.exp(2j * np.pi * rng.random((len(GAMMAS), c2.size)))
dt, T, ren = 2e-3, 200.0, 50
L1 = np.zeros(e1.shape); L4 = np.zeros(e1.shape)
cos_min = np.ones(c2.size); r_min = np.full(c2.size, np.inf); r_max = np.zeros(c2.size)
steps = int(T / dt)
for s in range(steps):
    k1 = rhs(c2, c3, e1, e4)
    k2 = rhs(*(x + 0.5 * dt * k for x, k in zip((c2, c3, e1, e4), k1)))
    k3 = rhs(*(x + 0.5 * dt * k for x, k in zip((c2, c3, e1, e4), k2)))
    k4 = rhs(*(x + dt * k for x, k in zip((c2, c3, e1, e4), k3)))
    c2, c3, e1, e4 = (x + dt / 6 * (a + 2 * b + 2 * c + d) for x, a, b, c, d in zip((c2, c3, e1, e4), k1, k2, k3, k4))
    if s % ren == 0:
        a1, a4 = np.abs(e1), np.abs(e4)
        L1 += 2 * np.log(a1); L4 += 2 * np.log(a4)                 # log I
        e1 /= a1; e4 /= a4
        rcur = abs(c3) ** 2 / abs(c2) ** 2
        r_min = np.minimum(r_min, rcur); r_max = np.maximum(r_max, rcur)
        cos_min = np.minimum(cos_min, np.cos(2 * (np.angle(c3) - np.angle(c2))))
a1, a4 = np.abs(e1), np.abs(e4)
L1 += 2 * np.log(a1); L4 += 2 * np.log(a4)
L1 = (L1 / T).reshape((len(GAMMAS),) + R.shape); L4 = (L4 / T).reshape((len(GAMMAS),) + R.shape)
LmaxG = np.maximum(L1, L4)                                           # units: M = 2 (rates scale with M)
Lmax = LmaxG[0]; L1 = L1[0]
Dm1, k1_, g1, _ = frozen(I2, I3, np.cos(TH))
Dm4, k4_, g4, _ = frozen(I3, I2, np.cos(TH))
gfro = np.maximum(g1, g4)
trapped_exact = Lmax < 0.05
trapped_frozen = gfro == 0
res["dimer_grid"] = {"points": int(R.size), "fraction_trapped_exact": float(trapped_exact.mean()),
                     "fraction_trapped_frozen": float(trapped_frozen.mean()),
                     "agreement_exact_vs_frozen": float((trapped_exact == trapped_frozen).mean()),
                     "equal_antiphase_Lmax": float(Lmax[np.argmin(abs(rr - 1)), np.argmin(abs(tt - np.pi))]),
                     "single_mode_check_rate_over_I2_frozen": float(2 * np.sqrt(3))}
# single hot mode exact check: I3 -> 0 row (r = 1/6) at any angle: rate of endpoint 1 vs 2 sqrt(3) I2
res["dimer_grid"]["trapped_fraction_damped"] = {str(g): float((LmaxG[i] < 0).mean()) for i, g in enumerate(GAMMAS[:, 0]) if g > 0}
res["dimer_grid"]["L1_r_min_over_frozen"] = float(np.median(L1[0] / g1[0]))
print(json.dumps(res, indent=1))
json.dump({**res, "grid_r": rr.tolist(), "grid_th": tt.tolist(), "Lmax": Lmax.tolist(), "gfrozen": gfro.tolist()},
          open(f"{O}/theory_trap.json", "w"))
json.dump({"LmaxG": LmaxG.tolist(), "gammas": GAMMAS[:, 0].tolist()}, open(f"{O}/theory_trap_damped.json", "w"))
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
cols = {0.02: "c", 0.05: "lime", 0.1: "w"}
for a, Z, title in [(ax[0], Lmax, "undamped: max endpoint Lyapunov exponent along the dimer orbit"),
                    (ax[1], gfro, "frozen-pair rate 2 sqrt(kappa^2 - Delta^2) (max of both ends)")]:
    m = a.pcolormesh(tt, np.log(rr), Z, shading="auto", cmap="magma")
    for i, g in enumerate(GAMMAS[:, 0]):
        if g > 0:
            a.contour(tt, np.log(rr), (LmaxG[i] < 0).astype(float), levels=[0.5], colors=cols[g], linewidths=1)
    plt.colorbar(m, ax=a)
    a.set_xlabel(r"initial $\theta_3=2(\varphi_3-\varphi_2)$"); a.set_ylabel(r"initial $\log(I_3/I_2)$")
    a.set_title(title + "\n(I2+I3 = 2; contours: damped exponent = 0 for gamma = 0.02 cyan, 0.05 green, 0.1 white)", fontsize=7)
fig.tight_layout(); fig.savefig(f"{O}/fig_theory_trap.png", dpi=150)
