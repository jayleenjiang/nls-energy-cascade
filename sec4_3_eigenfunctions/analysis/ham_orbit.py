# Deterministic (gamma = 0) Hamiltonian orbit from each ball center: distance to the start vs time.
import numpy as np, json
from scipy.integrate import solve_ivp
PI = np.pi
def wrap(a): return (a + PI) % (2 * PI) - PI
def rhs(t, u):
    p = u[:3] + 1j * u[3:]; I = np.abs(p) ** 2; M = I.sum()
    nb = np.zeros(3, complex); nb[0] = p[1] ** 2; nb[1] = p[0] ** 2 + p[2] ** 2; nb[2] = p[1] ** 2
    F = (2 * M - I) * p + 2 * nb * np.conj(p)          # same force as sampler (Re/Im parts)
    return np.concatenate([-F.imag, F.real])           # dx = -F_imag, dy = +F_real
def red(u):
    p = u[:3] + 1j * u[3:]; ph = np.angle(p)
    return np.array([*(np.abs(p) ** 2), wrap(2 * (ph[0] - ph[1])), wrap(2 * (ph[2] - ph[1]))])
def start(x0):
    ph = [x0[3] / 2, 0, x0[4] / 2]; p = np.sqrt(x0[:3]) * np.exp(1j * np.array(ph))
    return np.concatenate([p.real, p.imag])
C = {"c0": (1, 1, 1, 0, 0), "c1": (1, 1, 1, PI, PI), "c2": (2, 1, 1, 0, 0), "c3": (1, 1, 2, 0, 0),
     "c4": (1, 1, 1, PI / 2, -PI / 2), "c5": (1, 2, 1, PI, PI)}
out = {}
t = np.linspace(0, 5, 5001)
for k, x0 in C.items():
    x0 = np.array(x0, float)
    sol = solve_ivp(rhs, (0, 5), start(x0), t_eval=t, rtol=1e-10, atol=1e-12, method="DOP853")
    X = np.array([red(u) for u in sol.y.T])
    d = np.sqrt(((X[:, :3] - x0[:3]) ** 2).sum(1) + (wrap(X[:, 3:] - x0[3:]) ** 2).sum(1))
    mins = [i for i in range(1, len(t) - 1) if d[i] < d[i - 1] and d[i] < d[i + 1] and t[i] > 0.05]
    M = x0[:3].sum()
    out[k] = dict(t=t.tolist(), d=d.tolist(), X=X[::10].tolist())
    print(f"{k} M={M:.0f}: dtheta/dt at start = {np.round((X[1,3:]-X[0,3:])/(t[1]-t[0]),2)}; "
          f"closest returns (t, dist): {[(round(t[i],2), round(d[i],2)) for i in mins[:6]]}")
json.dump(out, open("ham_orbit.json", "w"))
