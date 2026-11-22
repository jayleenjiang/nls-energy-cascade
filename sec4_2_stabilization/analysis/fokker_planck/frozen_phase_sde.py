"""Monte Carlo cross-check of frozen_phase_fp.py: Euler--Maruyama for the two-angle diffusion,
L lanes, burn-in then time average of I1 sin th1 + I3 sin th3; SE over 16 lane groups."""
import sys, numpy as np
from frozen_phase_fp import drift, half_var
def run(model, I1, I3, gamma=0.1, T1=2.0, T3=8.0, dt=5e-4, burn=20.0, meas=40.0, L=16384, seed=0):
    rng = np.random.default_rng(seed)
    b1, b3 = drift(model, I1, I3); d1, d3 = half_var(model, I1, I3, gamma, T1, T3)
    s1, s3 = np.sqrt(2 * d1 * dt), np.sqrt(2 * d3 * dt)
    a = rng.uniform(-np.pi, np.pi, L); c = rng.uniform(-np.pi, np.pi, L)
    acc = np.zeros(L); n = 0
    for k in range(int((burn + meas) / dt)):
        da, dc = b1(a, c) * dt, b3(a, c) * dt
        a = a + da + s1 * rng.standard_normal(L); c = c + dc + s3 * rng.standard_normal(L)
        if k * dt >= burn and k % 10 == 0:
            acc += I1 * np.sin(a) + I3 * np.sin(c); n += 1
    g = (acc / n).reshape(16, -1).mean(1)
    return g.mean(), g.std(ddof=1) / 4
if __name__ == "__main__":
    for model, I1, I3 in [("canonical", 3.0, 3.0), ("hlns", 3.0, 3.0), ("canonical", 0.5, 2.0)]:
        for dt in [1e-3, 5e-4]:
            m, se = run(model, I1, I3, dt=dt)
            print(model, I1, I3, "dt", dt, "MC F =", round(m, 4), "+-", round(se, 4), flush=True)
